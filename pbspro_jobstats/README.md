One thing we don't might not bother with carefully on HPC is: 'how much memory (or CPUs) do I need to request?'

This little Python watcher that parses the stats and (optionally) makes plots of job resource usage over time (provided adequate Python3 libs installed--e.g. `pip3 install matplotlib`)

This is nicer to look at than `qstat -f | grep resource`

Most of us are probably requesting a lot of memory to be safe--but when is it too much? There are a lot of jobs requesting over 100 to 200+GB now on various HPC queues, leading to others' jobs having trouble finding resources. Necessary? Maybe....

This can work over SSH if you have SSH keys set up and your remote .bashrc $PATH set to include /usr/pbs/bin (see call examples under JobLogger:getstats())
