#mpi.sh
#!/bin/sh
#PBS−Nmpi #任务程序为mpi
#PBS−lnodes=4 # 需要四个节点来计算
pssh −h $PBS_NODEFILE mkdir −p /home/sTest 1>&2 # 在分配到的4个计算节点上创建对应路径
scp master:/home/sTest/mpi /home/sTest
#mpi为编译之后的可执行程序
pscp −h $PBS_NODEFILE /home/sTest/mpi /home/sTest 1>&2
#获取master节点中的mpi可执行程序，并分发到每个计算节点
mpiexec −np 4 −machinefile $PBS_NODEFILE /home/sTest/mpi
#在4个计算节点上运行可执行程序mpi