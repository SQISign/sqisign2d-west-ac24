add_executable(fp2.test_${SVARIANT_LOWER} ${ECX_DIR}/test/fp2-test.c)
	target_include_directories(fp2.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include ../include ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON})
	target_link_libraries(fp2.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}})
		
add_executable(poly-mul.test_${SVARIANT_LOWER} ${ECX_DIR}/test/poly-mul-test.c)
	target_include_directories(poly-mul.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON})
	target_link_libraries(poly-mul.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})
		
add_executable(poly-redc.test_${SVARIANT_LOWER} ${ECX_DIR}/test/poly-redc-test.c)
	target_include_directories(poly-redc.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include  ${INC_EC} ${INC_COMMON})
	target_link_libraries(poly-redc.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})


add_test(ec_fp2.test_${SVARIANT_LOWER} fp2.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_poly-mul.test_${SVARIANT_LOWER} poly-mul.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_poly-redc.test_${SVARIANT_LOWER} poly-redc.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
