 #!/bin/bash
for i in $(seq 1 1 20 ); do
rm output.sas
rm output
  if [ "$i" -lt "10" ];then
    ../../translate/translate.py ../benchmarks/openstacks-opt11-strips/p0$((i))-domain.pddl ../benchmarks/openstacks-opt11-strips/p0$((i)).pddl
    ../../preprocess/preprocess < output.sas
    mv output OPENSTACKS-$((i))
  else
    ../../translate/translate.py ../benchmarks/openstacks-opt11-strips/p$((i))-domain.pddl ../benchmarks/openstacks-opt11-strips/p$((i)).pddl
    ../../preprocess/preprocess < output.sas
    mv output OPENSTACKS-$((i))
  fi
done
