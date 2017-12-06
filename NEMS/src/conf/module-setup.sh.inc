__ms_shell=`ps -p$$ -ocmd= | cut -d'/' -f 2`

if [[ -d /lfs3 ]] ; then
    # We are on NOAA Jet
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /scratch3 ]] ; then
    # We are on NOAA Theia
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /apps/lmod/lmod/init/$__ms_shell
    fi
    module purge
elif [[ -d /gpfs/hps && -e /etc/SuSE-release ]] ; then
    # We are on NOAA Luna or Surge
    if ( ! eval module help > /dev/null 2>&1 ) ; then
	source /opt/modules/default/init/$__ms_shell
    fi
    module use /opt/cray/ari/modulefiles/
    module use /gpfs/hps/nco/ops/nwprod/lib/modulefiles/
    module purge
elif [[ -d /dcom && -d /hwrf ]] ; then
    # We are on NOAA Tide or Gyre
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        source /usrx/local/Modules/default/init/$__ms_shell
    fi
    module purge
elif [[ -d /glade ]] ; then
    # We are on NCAR Yellowstone
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        . /usr/share/Modules/init/$__ms_shell
    fi
    module purge
elif [[ -d /lustre && -d /ncrc ]] ; then
    # We are on GAEA. 
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        # We cannot simply load the module command.  The GAEA
        # /etc/profile modifies a number of module-related variables
        # before loading the module command.  Without those variables,
        # the module command fails.  Hence we actually have to source
        # /etc/profile here.
        source /etc/profile
    fi
    module purge
elif [[ -d /nobackupp8 ]] ; then
    # We are on Pleiades.
    if ( ! eval module help > /dev/null 2>&1 ) ; then
        . /usr/local/lib/init/global.profile
    fi
    module purge
else
    echo WARNING: UNKNOWN PLATFORM 1>&2
fi

unset __ms_shell
unset __ms_ksh_test
unset __ms_bash_test
unset $__ms_function_name
unset __ms_function_name
