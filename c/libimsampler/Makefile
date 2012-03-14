INCLUDES=-I../mt19937ar/ -I./ -I../common/
# CFLAGS=-std=c99 -D_GNU_SOURCE -fPIC -pthread -fexceptions \
# 	-D_FILE_OFFSET_BITS=64 -fvisibility=hidden $(INCLUDES)

CFLAGS=-std=c99 -pthread -fexceptions -fPIC $(INCLUDES)

IM_SOURCES=../common/sampling.c IM_unif.c IM_ego.c IMUtils.c IMOutData.c \
	SAWHeaps.c ../common/logspace.c ../common/GraphNodes.c \
	../common/objectives.c ../mt19937ar/mt19937ar.c ../common/utils.c
IM_MEX_SOURCES=GraphNodes_mex.c mex/IMRun_mxutils.c $(IM_SOURCES)

ctypes: libimsampler_ctypes

IMRun_Expert.mex:
	mex CFLAGS="$(CFLAGS)" -v -g $(INCLUDES) mex/IMRun_Expert.c $(IM_MEX_SOURCES) -outdir mex/
	cp mex/IMRun*.mex* ../matlab

IM_CTYPES_SOURCES=ctypes/IM_ctypes.c ../common/GraphNodes_ctypes.c $(IM_SOURCES)
IM_CTYPES_OBJS=$(addsuffix .o, $(basename $(IM_CTYPES_SOURCES)))

libimsampler_ctypes: clean_common $(IM_CTYPES_OBJS)
	@ # This command is for Linux
	$(CC) -shared -Wl,-soname,libimsampler_ctypes.so -o libimsampler_ctypes.so $(IM_CTYPES_OBJS)
	cp libimsampler_ctypes.so ../../python/
	
	@ # This command is for OS/X
	@ # $(CC) -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0 -exported_symbols_list globalsyms.sym -o ctypes/libimsampler_ctypes.dylib $(IM_CTYPES_OBJS)
	@ # cp ctypes/libimsampler_ctypes.dylib ../../python/

.PHONY: clean_common clean

clean_common:
	@echo cleaning common files
	rm -f ../common/*.o

clean:
	@echo cleaning up
	rm -f *.o mex/*.o mex/*.mex* ../matlab/*.mex* ctypes/*.o ../mt19937ar/*.o \
		../common/*.o \
		ctypes/libimsampler_ctypes.dylib \
		../../python/libimsampler_ctypes.dylib

