make -j  
gs=501; nt=2000
gs=501; nt=200
gs=501; nt=10
source=--disable-source-point 
source=
rm rcv.bin.bck
mv rcv.bin rcv.bin.bck
./build/mwd_kernel --nx $gs  --ny $gs --nz $gs --nt $nt --target-kernel 0 --mwd-type 0 --target-ts 0 --verify 1 $source 
#FROM HATEM
#d=96; ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 96 --ny 96 --nz 96 --nt 10 --mwd-type 1 --target-ts 2 --thread-group-size 6 --thx 1 --thy 2 --thz 3 --threads 36 --verify 1 #--disable-source-point 
