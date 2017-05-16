make -j  
#gs=512; nt=8;    verify=1; num_threads=12;  tgs=6
gs=512; nt=2101; verify=0; num_threads=18; tgs=6
gs=512; nt=2101; verify=0; num_threads=12; tgs=6
#gs=512; nt=10; verify=1; num_threads=12; tgs=6
#gs=96;  nt=2101;    verify=1; num_threads=12;  tgs=6
source=--disable-source-point 
source=
rm rcv.bin.bck
mv rcv.bin rcv.bin.bck
rm rcv-dia-*.bin
export $OMP_NUM_THREADS=$num_threads
#DIAMOND
cmd="./build/mwd_kernel --nx $gs  --ny $gs --nz $gs --nt $nt --mwd-type 1 --target-ts 2 --verify $verify $source  --npx 1 --npy 1 --npz 1 --thread-group-size $tgs --thx 1 --thy 2 --thz 3 --threads $num_threads --n-tests 1"
#SPATIAL
#cmd="./build/mwd_kernel --nx $gs  --ny $gs --nz $gs --nt $nt --target-kernel 0 --mwd-type 0 --target-ts 0 --verify 1 $source"
echo $cmd
#gdb --ex run --args \
$cmd

#FROM HATEM
#d=96; ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 96 --ny 96 --nz 96 --nt 10 --mwd-type 1 --target-ts 2 --thread-group-size 6 --thx 1 --thy 2 --thz 3 --threads 36 --verify 1 #--disable-source-point 
