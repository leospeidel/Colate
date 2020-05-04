project (htslib)

#add_library(htslibShared INTERFACE kstring kbitset klist kfunc khash khash_str2int knetfile kroundup kseq ksort hts hts_defs hts_endian hts_log hts_os sam vcf vcf_sweep vcfutils bgzf cram hfile regidx synced_bcf_reader tbx thread_pool)												
add_library(htslibShared SHARED hts hts_defs hts_endian hts_log hts_os kstring kbitset klist kfunc khash khash_str2int knetfile kroundup	kseq ksort sam vcf vcf_sweep vcfutils	bgzf cram hfile regidx synced_bcf_reader tbx thread_pool)
set_target_properties(htslibShared PROPERTIES LINKER_LANGUAGE C)

#add_library(htslibShared INTERFACE)
#target_sources(htslibShared INTERFACE hts hts_defs hts_endian hts_log hts_os kstring kbitset klist kfunc khash khash_str2int knetfile kroundup	kseq ksort sam vcf vcf_sweep vcfutils	bgzf cram hfile regidx synced_bcf_reader tbx thread_pool)												
#target_include_directories(htslibShared INTERFACE ${HTSLIB_DIR})

#add_library(htslibStatic STATIC kstring kbitset klist kfunc khash khash_str2int knetfile kroundup kseq ksort hts hts_defs hts_endian hts_log hts_os sam vcf vcf_sweep vcfutils bgzf cram hfile regidx synced_bcf_reader tbx thread_pool)												
#set_target_properties(htslibStatic PROPERTIES LINKER_LANGUAGE C)
