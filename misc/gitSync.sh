printf "git push on mac; pull on hpcc"
git push
ssh -t rl44w@ghpcc06.umassrc.org 'cd github/bpipes; git pull'
