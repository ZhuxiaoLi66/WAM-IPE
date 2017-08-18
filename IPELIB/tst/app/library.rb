module Library

  # REQUIRED METHODS (CALLED BY DRIVER)

  def lib_build(env,prepkit)
    cmd="cd #{prepkit[:srcdir]} && #{env.build.cmd}"
    Thread.exclusive { ext(cmd,{:msg=>"Build failed"}) }
    prepkit
  end

  def lib_build_post(env,buildkit)
    buildkit
  end

  def lib_build_prep(env)
    srcdir=File.join(env.build.ddts_root,"src")
    cmd="rsync -a #{File.join($DDTSHOME,"..","src")}/ #{srcdir}"
    Thread.exclusive { ext(cmd,{:msg=>"Error copying 'src' to '#{srcdir}'"}) }
    rundir=File.join(env.build.ddts_root,"run")
    cmd="rsync -a --no-recursive #{File.join($DDTSHOME,"..","run")}/* #{rundir}"
    Thread.exclusive { ext(cmd,{:msg=>"Error copying 'run' to '#{rundir}'"}) }
    bindir=File.join(env.build.ddts_root,"bin")
    FileUtils.mkdir_p(bindir)
    FileUtils.mkdir_p(File.join(env.build.ddts_root,"include"))
    FileUtils.mkdir_p(File.join(env.build.ddts_root,"lib"))
    {:bindir=>bindir,:rundir=>rundir,:srcdir=>srcdir}
  end

  def lib_data_trillian(env)
    link_data("/mnt/lustre/lus0/space/madden/IPE/test-suite-data")
  end

  def lib_data_zeus(env)
    link_data("/scratch1/portfolios/NCEPDEV/swpc/noscrub/Paul.A.Madden/IPE/test-suite-data")
  end

  def lib_data_theia(env)
    link_data("/scratch3/NCEPDEV/swpc/noscrub/Naomi.Maruyama/ipe/test-suite-data")
  end

  def lib_outfiles(env,path)
    restrs=['(.*/)(plasma\d\d)','(.*/)(fort\.20\d\d)']
    res=restrs.map { |e| Regexp.new(e) }
    outs=[]
    Find.find(path) do |e|
      res.each do |re|
        m=re.match(e)
        unless m.nil?
          outs << [m[1],m[2]]
          break
        end
      end
    end
    outs
  end

  def lib_queue_del_cmd(env)
    'qdel'
  end

  def lib_run(env,prepkit)
    rundir=prepkit
    datadir=valid_dir(File.join(tmp_dir,"data"))
    da="IPEDATA=#{datadir}"
    qu="IPEQUEUE=batch"
    wi="IPEWIND=#{datadir}"
    ma=env.run.machine
    co=env.run.compiler
    pa=env.run.parallelism
    ta=env.run.tasks
    cmd="cd #{rundir} && #{da} #{qu} #{wi} ./qsubipe #{ma} #{co} #{pa} #{ta}"
    logd "Submitting job with command: #{cmd}"
    output,status=Thread.exclusive do
      ext(cmd,{:msg=>"ERROR: Job submission failed"})
    end
    re1=Regexp.new('The job (\d+).* has been submitted.')
    re2=Regexp.new('Created run directory: (.*)')
    jobid=nil
    subdir=nil
    output.each do |e|
      e.chomp!
      logd e
      jobid=e.gsub(re1,'\1') if re1.match(e)
      subdir=e.gsub(re2,'\1') if re2.match(e)
    end
    die "ERROR: Job ID not found in submit-command output" if jobid.nil?
    job_activate(jobid,self)
    die "ERROR: Run directory not found in submit-command output" if subdir.nil?
    logi "Queued with job ID #{jobid}"
    outputdir=valid_dir(File.join(rundir,subdir))
    invoke(:wait_for_job,:run,env,jobid,outputdir)
    job_deactivate(jobid)
    outputdir
  end

  def lib_run_check(env,postkit)
    outputdir=postkit
    stdout=valid_file(File.join(outputdir,"output"))
    return outputdir if job_check(stdout,'IPE completed successfully')
    nil
  end

  def lib_run_post(env,runkit)
    outputdir=runkit
  end

  def lib_run_prep(env)
    uniq=env.run.ddts_root
    run=env.build.ddts_result[:rundir]
    logd "Copying #{run} -> #{uniq}"
    FileUtils.cp_r(run,uniq)
    bin=env.build.ddts_result[:bindir]
    logd "Linking #{bin} -> #{uniq}"
    FileUtils.ln_s(bin,uniq)
    rundir=File.join(uniq,File.basename(run))
    modcmd=valid_file(File.join(env.build.ddts_root,"src","modcmd"))
    logd "Copying #{modcmd} -> #{rundir}"
    FileUtils.cp(modcmd,rundir)
    nlfile=valid_file(File.join(rundir,'SMSnamelist'))
    mod_namelist_file(nlfile,env.run.namelists)
    rundir
  end

  # CUSTOM METHODS (NOT CALLED BY DRIVER)

  def link_data(dir)
    validate_data(dir)
    link=File.join(tmp_dir,"data")
    FileUtils.rm_f(link)
    FileUtils.ln_s(dir,link)
  end

  def mod_namelist_file(nlfile,nlenv)
    h=convert_o2h(nlenv)
    sets=h.reduce([]) do |m0,(n,kv)|
      inner=kv.reduce([]) do |m1,(k,v)|
        v="\"#{quote_string(v)}\""
        logd "Set namelist #{n}:#{k}=#{v}"
        m1.push("-s #{n}:#{k}=#{v}")
      end
      m0.concat(inner)
    end
    nml=valid_file(File.expand_path(File.join($DDTSHOME,"nml")))
    cmd="#{nml} -i #{nlfile} -o #{nlfile} #{sets.join(" ")}"
    Thread.exclusive { ext(cmd,{:msg=>"Failed to edit #{nlfile}"}) }
  end

  def validate_data(dir)
    logd "Validating data..."
    expected={
      'ipe_grid'   => '037327a5c67ca47e33c304966fe4ce12',
      'plasma00'   => '103d1974df32b5084d6697bd44e05b1d',
      'plasma01'   => '236305d3a437c7195392994a48ad4fe3',
      'plasma02'   => '818360354487d808e1285b3b364902e4',
      'plasma03'   => '7918b93f407d0f5d8300b2252d101145',
      'plasma04'   => '6ee788632b3d3d7773b7496b8fd2d535',
      'plasma05'   => '236182b31bbe0d019b2033e93f2834c5',
      'plasma06'   => 'd2a0b3ea7dfd092683cb203a4134ca54',
      'plasma07'   => '246ed6ad760085b6abc56b76ac549b52',
      'plasma08'   => 'ce9fed92dd54fc835e3edb6d807fe915',
      'plasma09'   => '1ab6bbbbdb77e28c7e70add4b6df2379',
      'plasma10'   => '070f4cbd83e896ecda6ca13c300e198f',
      'plasma11'   => '070f4cbd83e896ecda6ca13c300e198f',
      'ut_input'   => '282e3bfa655cd85329faeca13efeb6da',
      'ut_rec'     => '926b18e8268259bf96d0f1dcee7a4c31',
      'wind_input' => '5d88f51319c8f919c63e71b850ed91e9'
    }
    actual=Dir.glob("#{dir}/**/*")
    actual=actual.delete_if { |e| File.directory?(e) }
    actual=actual.reduce({}) do
      |m,e| m.merge!({e.sub(/#{dir}\/?/,'')=>Digest::MD5.file(e).to_s})
    end
    actual.keys.sort.each do |k|
      die "Unexpected data file: #{k}" unless expected[k]
      unless actual[k]==expected[k]
        logd "Checksum validation failed for #{dir}/#{k}"
        logd "  Expected #{expected[k]}"
        logd "    Actual #{actual[k]}"
        die "Error validating test-suite data, see #{logfile}"
      end
      logd "  #{k}: OK"
    end
    logd "Validating data: OK"
  end

  def wait_for_job_trillian(env,jobid,outputdir)
    ok=%w[E H M Q R S T W]
    begin
      sleep 30
      cmd="qstat -f #{jobid}"
      output,status=Thread.exclusive { ext(cmd,{:die=>false,:out=>false}) }
      live=false
      output.each do |e|
        logd e
        live=true if ok.include?(e.chomp.sub(/^ *job_state = (.)$/,'\1'))
      end
    end while live or not File.exist?(File.join(outputdir,"output.batch"))
  end

  def wait_for_job_zeus(env,jobid,outputdir)
    ok=%w[E H Q R S T W]
    # 'tolerance' is the number of seconds the batch system retains information
    # about completed jobs. If this interval passes without a non-error response
    # to queries, we may never receive confirmation that the job completed.
    # Consider this a batch-system failure and abort.
    tolerance=600
    batch_failure=false
    last_response=Time.now
    begin
      sleep 30
      tolog=[]
      cmd="qstat -f #{jobid}"
      output,status=Thread.exclusive { ext(cmd,{:die=>false,:out=>false}) }
      if status==0
        live=false
        last_response=Time.now
        output.each do |e|
          tolog.push(e)
          live=true if ok.include?(e.chomp.sub(/^ *job_state = (.)$/,'\1'))
        end
      else
        live=true
        now=Time.now
        logd "#{cmd} set error status #{status} at #{now}"
        if now-last_response > tolerance
          batch_failure=true
          live=false
        end
      end
    end while live
    logd "* Final batch info for job ID #{jobid}:"
    logd "--"
    tolog.each { |e| logd e }
    logd "--"
    die "* Batch system unresponsive for #{tolerance} seconds" if batch_failure
    re=Regexp.new('^ *exit_status = (\d+)')
    tolog.each do |e|
      m=re.match(e)
      return m[1].to_i if m
    end
    'unknown'
  end


  def wait_for_job_theia(env,jobid,outputdir)
    ok=%w[E H Q R S T W]
    # 'tolerance' is the number of seconds the batch system retains information
    # about completed jobs. If this interval passes without a non-error response
    # to queries, we may never receive confirmation that the job completed.
    # Consider this a batch-system failure and abort.
    tolerance=600
    batch_failure=false
    last_response=Time.now
    begin
      sleep 30
      tolog=[]
      cmd="qstat -f #{jobid}"
      output,status=Thread.exclusive { ext(cmd,{:die=>false,:out=>false}) }
      if status==0
        live=false
        last_response=Time.now
        output.each do |e|
          tolog.push(e)
          live=true if ok.include?(e.chomp.sub(/^ *job_state = (.)$/,'\1'))
        end
      else
        live=true
        now=Time.now
        logd "#{cmd} set error status #{status} at #{now}"
        if now-last_response > tolerance
          batch_failure=true
          live=false
        end
      end
    end while live
    logd "* Final batch info for job ID #{jobid}:"
    logd "--"
    tolog.each { |e| logd e }
    logd "--"
    die "* Batch system unresponsive for #{tolerance} seconds" if batch_failure
    re=Regexp.new('^ *exit_status = (\d+)')
    tolog.each do |e|
      m=re.match(e)
      return m[1].to_i if m
    end
    'unknown'
  end

end
