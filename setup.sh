cd python/cpp/
make clean
make depend
make

cd ../../
cd c/libsardsampler/
make

cd ../libimsampler/
make

cd ../../swendsen_wang/
mex sw_allall_flip_conditionals.c

cd ../python
python setup_config.py

cd ../
cp adaptive_saw.h5 ./python/