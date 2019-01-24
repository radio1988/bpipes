## TEST.sh
sleep 2 && echo two &
sleep 5 && echo five &
wait
echo First Round Finished At Five
sleep 5 && echo Ten &
wait
echo Second Round Finished At Ten
