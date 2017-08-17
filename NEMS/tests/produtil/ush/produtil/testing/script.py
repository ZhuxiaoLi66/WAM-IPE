"""!Takes an object tree from produtil.parse.Parser and turns it into
a flat bash script that will run the entire test suite, one test at a
time."""

import StringIO
 
__all__=['bash_functions','BashRunner']

import produtil.testing.parsetree

bash_functions=r'''
function report_start() {
  local parent
  parent=$( dirname "$1" )
  for x in 1 2 3 4 5 6 7 8 9 0 ; do
    if [[ ! -d "$parent" ]] ; then
      mkdir -p "$parent" || true
    fi
  done
  rt__TEST_REPORT_FILE="$1"
  rt__TEST_SUCCESS="YES"
  shift 1
  echo "$*" > "$rt__TEST_REPORT_FILE"
  date >> "$rt__TEST_REPORT_FILE"
}

function report_line() {
  echo "$*" >> "$rt__TEST_REPORT_FILE"
}

function report_stdin() {
  cat >> "$rt__TEST_REPORT_FILE"
}

function report_failure() {
  echo "$*" >> "$rt__TEST_REPORT_FILE"
  rt__TEST_SUCCESS=NO
}

function report_finish() {
  if [[ "$rt__TEST_SUCCESS" == "YES" ]] ; then
    echo "TEST PASSED AT $( date )" >> "$rt__TEST_REPORT_FILE"
  else
    echo "TEST FAILED AT $( date )" >> "$rt__TEST_REPORT_FILE"
    exit 1
  fi
}

function deliver_file() {
  local src tgt parent
  set -e
  src="$1"
  tgt="$2"
  if [[ -d "$tgt" ]] ; then
    tgt="$tgt"/$( basename "$src" )
  fi
  parent=$( dirname "$tgt" )
  if [[ ! -d "$parent" ]] ; then
    mkdir -p "$parent"
  fi
  tmpfile="$tgt.$$.$RANDOM.$RANDOM"
  cp -fp "$src" "$tmpfile"
  mv -fT "$tmpfile" "$tgt"
}

function bitcmp() {
  local src tgt bn result origtgt
  set -e
  src="$1"
  tgt="$2"
  if [[ -d "$tgt" ]] ; then
    bn=$( basename "$src" )
    tgt="$tgt/$bn"
  elif [[ -d "$src" ]] ; then
    bn=$( basename "$tgt" )
    src="$src/$bn"
  fi
  set +e
  if [[ ! -e "$src" ]] ; then
    report_failure "$tgt: MISSING BASELINE FILE"
    return 1
  elif [[ ! -e "$tgt" ]] ; then
    report_failure "$src: MISSING OUTPUT FILE"
    return 1
  fi
  cmp "$src" "$tgt"
  result=$?
  if [[ "$result" != 0 ]] ; then
    report_failure "$src: MISMATCH $tgt"
  else
    report_line "$src: bit-for-bit identical"
  fi
  return $result
}

function atparse {
    set +x # There is too much output if "set -x" is on.
    set +u # We expect empty variables in this function.
    set +e # We expect some evals to fail too.
    # Use __ in names to avoid clashing with variables in {var} blocks.
    local __text __before __after __during
    for __text in "$@" ; do
        if [[ $__text =~ ^([a-zA-Z][a-zA-Z0-9_]*)=(.*)$ ]] ; then
            eval "local ${BASH_REMATCH[1]}"
            eval "${BASH_REMATCH[1]}="'"${BASH_REMATCH[2]}"'
        else
            echo "ERROR: Ignoring invalid argument $__text\n" 1>&2
        fi
    done
    while IFS= read -r __text ; do
        while [[ "$__text" =~ ^([^@]*)(@\[[a-zA-Z_][a-zA-Z_0-9]*\]|@\[\'[^\']*\'\]|@\[@\]|@)(.*) ]] ; do
            __before="${BASH_REMATCH[1]}"
            __during="${BASH_REMATCH[2]}"
            __after="${BASH_REMATCH[3]}"
#            printf 'PARSE[%s|%s|%s]\n' "$__before" "$__during" "$__after"
            printf %s "$__before"
            if [[ "$__during" =~ ^@\[\'(.*)\'\]$ ]] ; then
                printf %s "${BASH_REMATCH[1]}"
            elif [[ "$__during" == '@[@]' ]] ; then
                printf @
            elif [[ "$__during" =~ ^@\[([a-zA-Z_][a-zA-Z_0-9]*)\] ]] ; then
                set -u
                eval 'printf %s "$'"${BASH_REMATCH[1]}"'"'
                set +u
            else
                printf '%s' "$__during"
            fi
            if [[ "$__after" == "$__text" ]] ; then
                break
            fi
            __text="$__after"
        done
        printf '%s\n' "$__text"
    done
}
'''

class MPICHRunner(produtil.testing.parsetree.Context):
    def __init__(self,scopes,token,run_mode,logger):
        super(MPICHRunner,self).__init__(
            scopes,token,run_mode,logger)
    def mpirunner(self,spawnProcess):
        out=StringIO.StringIO()
        out.write('mpirun')
        for rank in spawnProcess.iterrank():
            out.write(' -np %d %s'%(
                    rank.ranks(self),
                    ' '.join([r.bash_context(self)
                              for r in rank.args])))
        ret=out.getvalue()
        out.close()
        return ret

class LSFRunner(produtil.testing.parsetree.Context):
    def __init__(self,scopes,token,run_mode,logger):
        super(LSFRunner,self).__init__(
            scopes,token,run_mode,logger)
    def mpirunner(self,spawnProcess):
        prior=None
        for rank in spawnProcess.iterrank():
            prog=' '.join([r.bash_context(self)
                           for r in rank.args])
            if prior is not None and prior!=prog:
                raise NotImplementedError(self.error(
                        'MPMD is not yet supported for LSF.'))
        return 'mpirun.lsf '+prog

def runner_context_for(con):
    MPI=con.scopes[-1].resolve('plat%MPI').string_context(con)
    if MPI=='LSF':
        return LSFRunner(con.scopes,con.token,con.run_mode,con.logger)
    elif MPI=='MPICH':
        return MPICHRunner(con.scopes,con.token,con.run_mode,con.logger)
    else:
        raise NotImplementedError(con.error(
                'Unknown or unsupported MPI implementation "%s"'%(
                    elipses(MPI),)))

class BashRunner(object):
    def __init__(self):
        super(BashRunner,self).__init__()
    def make_runner(self,parser,output_file,dry_run=False,
                    setarith=None):
        logger=parser.logger
        runset=parser.setarith(setarith)
        logger.info('%s: generate bash script'%(output_file,))
        if dry_run: return
        with open(output_file,'wt') as out:
            out.write(r'''#! /usr/bin/env bash

%s

set -xe

'''%(bash_functions,))

            seen=False
            for runcon in runset:
                runme,con=runcon.as_tuple
                seen=True
                out.write(runme.bash_context(runner_context_for(con)))
                out.write("\n\n")
        if not seen:
            raise ValueError('ERROR: No "run" statments seen; nothing to do.\n');
