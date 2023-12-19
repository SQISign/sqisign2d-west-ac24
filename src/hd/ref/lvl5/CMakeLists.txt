set(SOURCE_FILES_HD_${SVARIANT_UPPER}_REF
    # ${ECX_DIR}/poly-mul.c 
    # ${ECX_DIR}/poly-redc.c 
    # ${ECX_DIR}/ec.c 
    # ${ECX_DIR}/tedwards.c 
    # ${ECX_DIR}/kps.c 
    # ${ECX_DIR}/xisog.c 
    # ${ECX_DIR}/xeval.c 
    # ${ECX_DIR}/isog_chains.c 
    # ${ECX_DIR}/basis.c
    ${HDX_DIR}/hd.c 
)

add_library(${LIB_HD_${SVARIANT_UPPER}} ${SOURCE_FILES_HD_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_HD_${SVARIANT_UPPER}} PRIVATE ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_EC})
target_compile_options(${LIB_HD_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
