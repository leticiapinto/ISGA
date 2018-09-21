for k in 1 2 ; do
	for j in 1 2 3 4 5 6 7 8 9 10 ; do
		result=$((10*${k} + ${j} -10))
		echo ${result}
		for i in 1 2 3 4 5 6 7 8 9 10 ; do
			./ig.out  ../Datasets/1000_25/1000_25_${i}.dat >> 1000_25_${i}.txt
			./ig.out  ../Datasets/1000_50/1000_50_${i}.dat >> 1000_50_${i}.txt
			./ig.out  ../Datasets/1000_75/1000_75_${i}.dat >> 1000_75_${i}.txt
			./ig.out  ../Datasets/1000_100/1000_100_${i}.dat >> 1000_100_${i}.txt

			./ig.out  ../Datasets/2000_25/2000_25_${i}.dat >> 2000_25_${i}.txt
			./ig.out  ../Datasets/2000_50/2000_50_${i}.dat >> 2000_50_${i}.txt
			./ig.out  ../Datasets/2000_75/2000_75_${i}.dat >> 2000_75_${i}.txt
			./ig.out  ../Datasets/2000_100/2000_100_${i}.dat >> 2000_100_${i}.txt
		done
	done
done