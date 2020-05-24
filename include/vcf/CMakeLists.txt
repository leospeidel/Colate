project (htslib)

#set (CMAKE_SHARED_LINKER_FLAGS "-lz lm lbz2 -lzma -lcurl -lssl,--as-needed")
add_library(hts SHARED
bcf_sr_sort.h          header.h               hts_internal.h         textutils_internal.h   version.h
config.h               hfile_internal.h       sam_internal.h         thread_pool_internal.h

bcf_sr_sort.c          errmod.c               hfile_gcs.c            hfile_s3_write.c    kfunc.c             
multipart.c            regidx.c               vcf.c                  bgzf.c              faidx.c             
hfile_libcurl.c        hts.c                  knetfile.c             region.c            tbx.c               
vcf_sweep.c            header.c               hfile_net.c            hts_os.c            kstring.c           
probaln.c              sam.c                  textutils.c            vcfutils.c          hfile.c             
hfile_s3.c             md5.c                  realn.c                synced_bcf_reader.c thread_pool.c

cram/cram.h            cram/cram_encode.h     cram/cram_samtools.h   cram/mFILE.h           cram/os.h              cram/rANS_static.h
cram/cram_codecs.h     cram/cram_index.h      cram/cram_stats.h      cram/misc.h            cram/pooled_alloc.h    cram/string_alloc.h
cram/cram_decode.h     cram/cram_io.h         cram/cram_structs.h    cram/open_trace_file.h cram/rANS_byte.h	
cram/cram_codecs.c     cram/cram_encode.c     cram/cram_index.c      cram/cram_samtools.c   cram/mFILE.c           cram/pooled_alloc.c    cram/string_alloc.c
cram/cram_decode.c     cram/cram_external.c   cram/cram_io.c         cram/cram_stats.c      cram/open_trace_file.c cram/rANS_static.c

htslib/bgzf.h              htslib/hts.h               htslib/hts_os.h            htslib/khash_str2int.h     htslib/ksort.h             htslib/synced_bcf_reader.h htslib/vcf_sweep.h
htslib/cram.h              htslib/hts_defs.h          htslib/kbitset.h           htslib/klist.h             htslib/kstring.h           htslib/tbx.h               htslib/vcfutils.h
htslib/faidx.h             htslib/hts_endian.h        htslib/kfunc.h             htslib/knetfile.h          htslib/regidx.h            htslib/thread_pool.h
htslib/hfile.h             htslib/hts_log.h           htslib/khash.h             htslib/kseq.h              htslib/sam.h               htslib/vcf.h

os/lzma_stub.h
htslib.hpp htslib.cpp
)
target_link_libraries(hts OpenSSL::Crypto z m bz2 lzma curl)

