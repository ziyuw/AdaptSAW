INCLUDES=-I../mt19937ar/ -I./ -I../common/
#CFLAGS=-std=c99 -D_GNU_SOURCE -fPIC -pthread -m32 -fexceptions -D_FILE_OFFSET_BITS=64 -fvisibility=hidden $(INCLUDES)
CFLAGS=-std=c99 -pthread -fexceptions -fPIC $(INCLUDES)

SARD_SOURCES=SARD_unif.c SARD_ego.c SARDOutData.c SARDUtils.c SAWHeaps.c \
	../mt19937ar/mt19937ar.c ../common/utils.c ../common/GraphNodes.c \
	../common/logspace.c ../common/sampling.c ../common/objectives.c

ctypes: libsardsampler_ctypes

SARD_CTYPES_SOURCES=ctypes/SARD_ctypes.c ../common/GraphNodes_ctypes.c $(SARD_SOURCES)
SARD_CTYPES_OBJS=$(addsuffix .o, $(basename $(SARD_CTYPES_SOURCES)))

libsardsampler_ctypes: clean_common $(SARD_CTYPES_OBJS)
	@ # This command is for Linux
	#$(CC) -dynamiclib -Wl,-soname,libsardsampler_ctypes.so -export-symbols globalsyms.sym -o libsardsampler_ctypes.so $(SARD_CTYPES_OBJS)
	$(CC) -shared -Wl,-soname,libsardsampler_ctypes.so -o libsardsampler_ctypes.so $(SARD_CTYPES_OBJS)	
	cp libsardsampler_ctypes.so ../../python/

	@ # This command is for OS/X
	@ # $(CC) -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0 -exported_symbols_list globalsyms.sym -o ctypes/libsardsampler_ctypes.dylib $(SARD_CTYPES_OBJS)
	@ # cp ctypes/libsardsampler_ctypes.dylib ../../python/

.PHONY: clean_common clean

clean_common:
	@echo cleaning common files
	rm -f ../common/*.o

clean:
	@echo cleaning up
	rm -f *.o ctypes/*.o ../mt19937ar/*.o ../common/*.o \
		  ctypes/libsardsampler_ctypes.dylib \
		  ../../python/libsardsampler_ctypes.dylib

