
set(SOURCE_FILES_ID2ISO_GENERIC_REF
    ${ID2ISOX_DIR}/id2iso.c
)

add_library(${LIB_ID2ISO_${SVARIANT_UPPER}} ${SOURCE_FILES_ID2ISO_GENERIC_REF})
target_include_directories(${LIB_ID2ISO_${SVARIANT_UPPER}} PRIVATE common ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_HD} ${INC_ID2ISO} ${INC_COMMON} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}})
target_link_libraries(${LIB_ID2ISO_${SVARIANT_UPPER}} ${LIB_KLPT_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})
target_compile_options(${LIB_ID2ISO_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
