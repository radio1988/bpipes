printf ">>>git push on mac; pull on hpcc\n"

printf "\n\n\n>>>Pushing\n"
git push

printf "\n\n\n>>>Pulling on HPCC\n"
ssh -t rl44w@ghpcc06.umassrc.org 'cd github/bpipes; git pull'
