add_executable(fp2.test_${SVARIANT_LOWER} ${ECX_DIR}/test/fp2-test.c)
	target_include_directories(fp2.test_${SVARIANT_LOWER} PUBLIC ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include ../include ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON})
	target_link_libraries(fp2.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}})
		
add_executable(poly-mul.test_${SVARIANT_LOWER} ${ECX_DIR}/test/poly-mul-test.c)
	target_include_directories(poly-mul.test_${SVARIANT_LOWER} PUBLIC ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON})
	target_link_libraries(poly-mul.test_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}})
		
add_executable(poly-redc.test_${SVARIANT_LOWER} ${ECX_DIR}/test/poly-redc-test.c)
	target_include_directories(poly-redc.test_${SVARIANT_LOWER} PUBLIC ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include  ${INC_EC} ${INC_COMMON})
	target_link_libraries(poly-redc.test_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}})

add_executable(basis-gen.test_${SVARIANT_LOWER} ${ECX_DIR}/test/basis-gen-test.c)
	target_include_directories(basis-gen.test_${SVARIANT_LOWER} PUBLIC ${ECX_DIR}/test ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} .)
	target_link_libraries(basis-gen.test_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_EC_${SVARIANT_UPPER}})

add_test(ec_fp2.test_${SVARIANT_LOWER} fp2.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_poly-mul.test_${SVARIANT_LOWER} poly-mul.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_poly-redc.test_${SVARIANT_LOWER} poly-redc.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_basis_gen.test_${SVARIANT_LOWER} basis-gen.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
