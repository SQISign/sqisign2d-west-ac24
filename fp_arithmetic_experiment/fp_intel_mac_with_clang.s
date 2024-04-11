	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 14, 0	sdk_version 14, 2
	.globl	_fp_add                         ## -- Begin function fp_add
	.p2align	4, 0x90
_fp_add:                                ## @fp_add
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	(%rdx), %r8
	movq	8(%rdx), %rax
	addq	(%rsi), %r8
	movq	16(%rdx), %rcx
	adcq	8(%rsi), %rax
	movq	24(%rdx), %rdx
	adcq	16(%rsi), %rcx
	adcq	24(%rsi), %rdx
	movq	%rdx, %rsi
	shrq	$59, %rsi
	movq	%rsi, %r9
	negq	%r9
	movabsq	$-360287970189639680, %r10      ## imm = 0xFB00000000000000
	andq	%r10, %r9
	addq	%r8, %rsi
	adcq	$0, %rax
	adcq	$0, %rcx
	adcq	%rdx, %r9
	movq	%r9, %rdx
	shrq	$59, %rdx
	movq	%rdx, %r8
	negq	%r8
	andq	%r10, %r8
	addq	%rsi, %rdx
	adcq	$0, %rax
	adcq	$0, %rcx
	adcq	%r9, %r8
	movq	%rdx, (%rdi)
	movq	%rax, 8(%rdi)
	movq	%rcx, 16(%rdi)
	movq	%r8, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_sub                         ## -- Begin function fp_sub
	.p2align	4, 0x90
_fp_sub:                                ## @fp_sub
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	xorl	%r9d, %r9d
	movq	(%rsi), %rcx
	subq	(%rdx), %rcx
	movq	8(%rsi), %rax
	movl	$0, %r8d
	sbbq	%r8, %r8
	sarq	$63, %r8
	subq	8(%rdx), %rax
	movl	$0, %r10d
	sbbq	%r10, %r10
	addq	%r8, %rax
	adcq	%r8, %r10
	sarq	$63, %r10
	movq	16(%rsi), %r8
	subq	16(%rdx), %r8
	movl	$0, %r11d
	sbbq	%r11, %r11
	addq	%r10, %r8
	adcq	%r10, %r11
	sarq	$63, %r11
	movq	24(%rsi), %rsi
	subq	24(%rdx), %rsi
	movl	$0, %edx
	sbbq	%rdx, %rdx
	addq	%r11, %rsi
	adcq	%r11, %rdx
	sarq	$63, %rdx
	movl	%edx, %r10d
	andl	$2, %r10d
	subq	%r10, %rcx
	sbbq	%r9, %r9
	sarq	$63, %r9
	addq	%r9, %rax
	adcq	$0, %r9
	sarq	$63, %r9
	addq	%r9, %r8
	adcq	$0, %r9
	shrq	$63, %r9
	movabsq	$-720575940379279360, %r10      ## imm = 0xF600000000000000
	andq	%rdx, %r10
	orq	%r9, %r10
	subq	%r10, %rsi
	movq	%rsi, %rdx
	shrq	$59, %rdx
	movq	%rdx, %r9
	negq	%r9
	movabsq	$-360287970189639680, %r10      ## imm = 0xFB00000000000000
	andq	%r9, %r10
	addq	%rcx, %rdx
	adcq	$0, %rax
	adcq	$0, %r8
	adcq	%rsi, %r10
	movq	%rdx, (%rdi)
	movq	%rax, 8(%rdi)
	movq	%r8, 16(%rdi)
	movq	%r10, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_neg                         ## -- Begin function fp_neg
	.p2align	4, 0x90
_fp_neg:                                ## @fp_neg
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	xorl	%r8d, %r8d
	movq	$-2, %rdx
	subq	(%rsi), %rdx
	sbbq	%r8, %r8
	sarq	$63, %r8
	movq	8(%rsi), %rax
	movq	16(%rsi), %rcx
	notq	%rax
	addq	%r8, %rax
	adcq	$0, %r8
	sarq	$63, %r8
	notq	%rcx
	addq	%r8, %rcx
	adcq	$0, %r8
	sarq	$63, %r8
	subq	24(%rsi), %r8
	movabsq	$720575940379279359, %rsi       ## imm = 0x9FFFFFFFFFFFFFF
	addq	%r8, %rsi
	movq	%rsi, %r8
	shrq	$59, %r8
	movq	%r8, %r9
	negq	%r9
	movabsq	$-360287970189639680, %r10      ## imm = 0xFB00000000000000
	andq	%r9, %r10
	addq	%rdx, %r8
	adcq	$0, %rax
	adcq	$0, %rcx
	adcq	%rsi, %r10
	movq	%r8, (%rdi)
	movq	%rax, 8(%rdi)
	movq	%rcx, 16(%rdi)
	movq	%r10, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_mul                         ## -- Begin function fp_mul
	.p2align	4, 0x90
_fp_mul:                                ## @fp_mul
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%rdi, -176(%rbp)                ## 8-byte Spill
	movq	(%rsi), %r9
	movq	8(%rsi), %r15
	movq	16(%rsi), %r10
	movq	24(%rsi), %rax
	movq	(%rdx), %r11
	movq	8(%rdx), %rdi
	movq	%rdi, -104(%rbp)                ## 8-byte Spill
	movq	16(%rdx), %rbx
	movq	24(%rdx), %rsi
	movq	%r11, %rdx
	mulxq	%r9, %rcx, %r8
	movq	%rcx, -120(%rbp)                ## 8-byte Spill
	movq	%rdi, %rdx
	mulxq	%r15, %rcx, %rdx
	movq	%rdx, -80(%rbp)                 ## 8-byte Spill
	movq	%r15, -96(%rbp)                 ## 8-byte Spill
	movq	%rbx, %rdx
	movq	%r10, -112(%rbp)                ## 8-byte Spill
	mulxq	%r10, %r14, %rdx
	movq	%r14, -56(%rbp)                 ## 8-byte Spill
	movq	%rdx, -64(%rbp)                 ## 8-byte Spill
	movq	%rsi, %rdx
	mulxq	%rax, %r14, %rdx
	movq	%r14, -88(%rbp)                 ## 8-byte Spill
	movq	%rdx, -72(%rbp)                 ## 8-byte Spill
	movq	%rax, %r14
	movq	%rax, -168(%rbp)                ## 8-byte Spill
	movq	%rdi, %rdx
	mulxq	%r9, %rax, %r12
	movq	%rsi, %rdx
	mulxq	%r9, %r13, %rdi
	addq	%r8, %rax
	adcq	%rcx, %r12
	adcq	-80(%rbp), %r13                 ## 8-byte Folded Reload
	movq	%rsi, %rdx
	mulxq	%r10, %r8, %rcx
	adcq	-56(%rbp), %rdi                 ## 8-byte Folded Reload
	adcq	-64(%rbp), %r8                  ## 8-byte Folded Reload
	adcq	-88(%rbp), %rcx                 ## 8-byte Folded Reload
	movq	%r11, %rdx
	mulxq	%r15, %rdx, %r10
	setb	-56(%rbp)                       ## 1-byte Folded Spill
	addq	%rax, %rdx
	movq	%rdx, -64(%rbp)                 ## 8-byte Spill
	adcq	%r12, %r10
	movq	%r11, %rdx
	mulxq	%r14, %rax, %r12
	adcq	%r13, %rax
	adcq	%rdi, %r12
	movq	%rbx, %rdx
	mulxq	%r14, %r15, %rdx
	adcq	%r8, %r15
	adcq	%rcx, %rdx
	movq	%rdx, -88(%rbp)                 ## 8-byte Spill
	movzbl	-56(%rbp), %ecx                 ## 1-byte Folded Reload
	adcq	-72(%rbp), %rcx                 ## 8-byte Folded Reload
	movq	%rcx, -144(%rbp)                ## 8-byte Spill
	movq	%rbx, %rdx
	mulxq	%r9, %rcx, %rdi
	addq	%r10, %rcx
	adcq	%rax, %rdi
	movq	%rsi, %rdx
	movq	-96(%rbp), %r8                  ## 8-byte Reload
	mulxq	%r8, %r14, %r13
	adcq	%r12, %r14
	adcq	%r15, %r13
	movq	%r11, %rdx
	movq	-112(%rbp), %r15                ## 8-byte Reload
	mulxq	%r15, %r10, %r9
	setb	-41(%rbp)                       ## 1-byte Folded Spill
	addq	%rcx, %r10
	movq	%r10, -72(%rbp)                 ## 8-byte Spill
	movq	%rbx, %rdx
	mulxq	%r8, %rsi, %rcx
	movq	%rsi, -128(%rbp)                ## 8-byte Spill
	movq	%rcx, -152(%rbp)                ## 8-byte Spill
	movq	-104(%rbp), %r11                ## 8-byte Reload
	movq	%r11, %rdx
	mulxq	%r15, %rax, %r12
	movq	%rax, -136(%rbp)                ## 8-byte Spill
	adcq	%rdi, %r9
	movq	%r9, -112(%rbp)                 ## 8-byte Spill
	movabsq	$360287970189639680, %rdi       ## imm = 0x500000000000000
	movq	-120(%rbp), %rcx                ## 8-byte Reload
	movq	%rcx, %rdx
	mulxq	%rdi, %rdx, %r8
	movq	%rdx, -80(%rbp)                 ## 8-byte Spill
	setb	%bl
	leaq	(%rsi,%rax), %r15
	addq	%r9, %r15
	movq	-64(%rbp), %rdx                 ## 8-byte Reload
	mulxq	%rdi, %rax, %r9
	leaq	(%rcx,%rcx,4), %rdx
	shlq	$56, %rdx
	leaq	(%rdx,%r15), %rcx
	movq	%rcx, -160(%rbp)                ## 8-byte Spill
	addq	%r8, %rax
	movq	%rax, -56(%rbp)                 ## 8-byte Spill
	movq	%r10, %rdx
	mulxq	%rdi, %rax, %rsi
	adcq	%r9, %rax
	movq	%rax, -96(%rbp)                 ## 8-byte Spill
	movq	%rcx, %rdx
	mulxq	%rdi, %r15, %rcx
	adcq	%rsi, %r15
	adcq	-144(%rbp), %rcx                ## 8-byte Folded Reload
	addb	$255, %bl
	movq	%r11, %rdx
	mulxq	-168(%rbp), %rbx, %rax          ## 8-byte Folded Reload
	adcq	%r14, %rbx
	adcq	%r13, %rax
	movq	%rax, -104(%rbp)                ## 8-byte Spill
	movzbl	-41(%rbp), %edi                 ## 1-byte Folded Reload
	adcq	-88(%rbp), %rdi                 ## 8-byte Folded Reload
	adcq	$0, %rcx
	movq	-128(%rbp), %rdx                ## 8-byte Reload
	addq	-136(%rbp), %rdx                ## 8-byte Folded Reload
	adcq	-152(%rbp), %r12                ## 8-byte Folded Reload
	setb	%r11b
	xorl	%r8d, %r8d
	movq	-120(%rbp), %rax                ## 8-byte Reload
	negq	%rax
	sbbq	%r8, %r8
	sarq	$63, %r8
	movq	%r8, %rsi
	subq	-64(%rbp), %rsi                 ## 8-byte Folded Reload
	sbbq	$0, %r8
	sarq	$63, %r8
	movq	%r8, %r9
	subq	-72(%rbp), %r9                  ## 8-byte Folded Reload
	sbbq	$0, %r8
	sarq	$63, %r8
	movq	-80(%rbp), %r14                 ## 8-byte Reload
	addq	%r8, %r14
	adcq	$0, %r8
	subq	-160(%rbp), %r14                ## 8-byte Folded Reload
	sbbq	$0, %r8
	sarq	$63, %r8
	movq	-56(%rbp), %r13                 ## 8-byte Reload
	addq	%r8, %r13
	adcq	$0, %r8
	sarq	$63, %r8
	movq	-96(%rbp), %r10                 ## 8-byte Reload
	addq	%r8, %r10
	adcq	$0, %r8
	sarq	$63, %r8
	addq	%r8, %r15
	adcq	$0, %r8
	sarq	$63, %r8
	addq	-112(%rbp), %rdx                ## 8-byte Folded Reload
	adcq	%rbx, %r12
	movzbl	%r11b, %r11d
	adcq	-104(%rbp), %r11                ## 8-byte Folded Reload
	adcq	$0, %rdi
	adcq	%rcx, %r8
	addq	-120(%rbp), %rax                ## 8-byte Folded Reload
	movq	-64(%rbp), %rcx                 ## 8-byte Reload
	adcq	$0, %rcx
	setb	%al
	addq	%rsi, %rcx
	movzbl	%al, %eax
	adcq	-72(%rbp), %rax                 ## 8-byte Folded Reload
	setb	%cl
	addq	%r9, %rax
	movzbl	%cl, %esi
	adcq	%rdx, %rsi
	setb	%dl
	xorl	%eax, %eax
	addq	%r13, %r12
	setb	%al
	xorl	%ecx, %ecx
	addq	%r10, %r11
	setb	%cl
	xorl	%r9d, %r9d
	addq	%r15, %rdi
	setb	%r9b
	addq	%r14, %rsi
	movzbl	%dl, %esi
	adcq	%rsi, %r12
	adcq	%r11, %rax
	adcq	%rdi, %rcx
	adcq	%r8, %r9
	movq	-176(%rbp), %rdx                ## 8-byte Reload
	movq	%r12, (%rdx)
	movq	%rax, 8(%rdx)
	movq	%rcx, 16(%rdx)
	movq	%r9, 24(%rdx)
	addq	$8, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_sqr                         ## -- Begin function fp_sqr
	.p2align	4, 0x90
_fp_sqr:                                ## @fp_sqr
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%rdi, -144(%rbp)                ## 8-byte Spill
	movq	(%rsi), %rax
	movq	8(%rsi), %r10
	movq	24(%rsi), %r8
	movq	%r10, %rdx
	mulxq	%rax, %r15, %r11
	movq	%r8, %rdx
	mulxq	%rax, %rbx, %rdi
	movq	16(%rsi), %r9
	movq	%r9, %rdx
	mulxq	%rax, %r12, %rcx
	addq	%r11, %r12
	adcq	%rbx, %rcx
	mulxq	%r10, %rbx, %rsi
	setb	%r14b
	leaq	(%rcx,%rbx), %rdx
	shldq	$1, %r12, %rdx
	movq	%rdx, %r11
	shldq	$1, %r15, %r12
	addq	%r15, %r15
	movq	%rax, %rdx
	mulxq	%rax, %rdx, %rax
	movq	%rdx, -80(%rbp)                 ## 8-byte Spill
	addq	%r15, %rax
	movq	%rax, -72(%rbp)                 ## 8-byte Spill
	movq	%r10, %rdx
	mulxq	%r10, %rdx, %rax
	adcq	%r12, %rdx
	movq	%rdx, -64(%rbp)                 ## 8-byte Spill
	setb	%dl
	adcq	%rax, %r11
	movq	%r11, -56(%rbp)                 ## 8-byte Spill
	setb	-42(%rbp)                       ## 1-byte Folded Spill
	addb	$255, %dl
	adcq	$0, %rax
	movq	%r8, %rdx
	mulxq	%r9, %r13, %r15
	setb	-41(%rbp)                       ## 1-byte Folded Spill
	addb	$255, %r14b
	mulxq	%r10, %r12, %rax
	adcq	%rdi, %r12
	adcq	%r13, %rax
	setb	%dl
	addq	%rbx, %rcx
	adcq	%rsi, %r12
	movzbl	%dl, %esi
	adcq	$0, %rax
	adcq	%r15, %rsi
	movq	%rsi, %r10
	shldq	$1, %rax, %rsi
	shldq	$1, %r12, %rax
	movq	%r9, %rdx
	mulxq	%r9, %rdx, %rdi
	shldq	$1, %rcx, %r12
	shrq	$63, %r10
	movq	%r10, -112(%rbp)                ## 8-byte Spill
	xorl	%r10d, %r10d
	addq	%rdx, %r12
	setb	%r10b
	xorl	%r11d, %r11d
	addq	%rax, %rdi
	movq	%rdi, -136(%rbp)                ## 8-byte Spill
	movq	%r8, %rdx
	mulxq	%r8, %rcx, %rax
	movq	%rax, -104(%rbp)                ## 8-byte Spill
	setb	%r11b
	xorl	%eax, %eax
	addq	%rsi, %rcx
	movq	%rcx, -128(%rbp)                ## 8-byte Spill
	setb	%al
	movq	%rax, -120(%rbp)                ## 8-byte Spill
	movq	-80(%rbp), %r15                 ## 8-byte Reload
	leaq	(%r15,%r15,4), %rax
	shlq	$56, %rax
	movq	-56(%rbp), %rcx                 ## 8-byte Reload
	leaq	(%rax,%rcx), %r8
	movabsq	$360287970189639680, %rax       ## imm = 0x500000000000000
	movq	%r15, %rdx
	mulxq	%rax, %r13, %rcx
	movq	%rcx, -96(%rbp)                 ## 8-byte Spill
	movq	-72(%rbp), %rdx                 ## 8-byte Reload
	mulxq	%rax, %r14, %rcx
	movq	%rcx, -88(%rbp)                 ## 8-byte Spill
	xorl	%ecx, %ecx
	negq	%r15
	sbbq	%rcx, %rcx
	sarq	$63, %rcx
	movq	%rcx, %r9
	subq	%rdx, %r9
	sbbq	$0, %rcx
	sarq	$63, %rcx
	movq	%rcx, %rdi
	movq	-64(%rbp), %rdx                 ## 8-byte Reload
	subq	%rdx, %rdi
	sbbq	$0, %rcx
	sarq	$63, %rcx
	addq	%rcx, %r13
	adcq	$0, %rcx
	subq	%r8, %r13
	sbbq	$0, %rcx
	sarq	$63, %rcx
	mulxq	%rax, %rsi, %rbx
	addq	-96(%rbp), %r14                 ## 8-byte Folded Reload
	adcq	-88(%rbp), %rsi                 ## 8-byte Folded Reload
	movq	%r8, %rdx
	mulxq	%rax, %rdx, %r8
	adcq	%rbx, %rdx
	adcq	-104(%rbp), %r8                 ## 8-byte Folded Reload
	addq	%rcx, %r14
	adcq	$0, %rcx
	sarq	$63, %rcx
	addq	%rcx, %rsi
	adcq	$0, %rcx
	sarq	$63, %rcx
	addq	%rcx, %rdx
	adcq	$0, %rcx
	addq	-112(%rbp), %r8                 ## 8-byte Folded Reload
	sarq	$63, %rcx
	addq	%rcx, %r8
	addb	$255, -42(%rbp)                 ## 1-byte Folded Spill
	movzbl	-41(%rbp), %ecx                 ## 1-byte Folded Reload
	adcq	%r12, %rcx
	adcq	-136(%rbp), %r10                ## 8-byte Folded Reload
	adcq	-128(%rbp), %r11                ## 8-byte Folded Reload
	movq	-120(%rbp), %rbx                ## 8-byte Reload
	adcq	%r8, %rbx
	addq	-80(%rbp), %r15                 ## 8-byte Folded Reload
	adcq	-72(%rbp), %r9                  ## 8-byte Folded Reload
	adcq	-64(%rbp), %rdi                 ## 8-byte Folded Reload
	movq	-56(%rbp), %r9                  ## 8-byte Reload
	adcq	$0, %r9
	setb	%r8b
	xorl	%eax, %eax
	addq	%r14, %rcx
	setb	%al
	xorl	%edi, %edi
	addq	%rsi, %r10
	setb	%dil
	xorl	%esi, %esi
	addq	%rdx, %r11
	setb	%sil
	addq	%r13, %r9
	movzbl	%r8b, %edx
	adcq	%rdx, %rcx
	adcq	%r10, %rax
	adcq	%r11, %rdi
	adcq	%rbx, %rsi
	movq	-144(%rbp), %rdx                ## 8-byte Reload
	movq	%rcx, (%rdx)
	movq	%rax, 8(%rdx)
	movq	%rdi, 16(%rdx)
	movq	%rsi, 24(%rdx)
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_n_sqr                       ## -- Begin function fp_n_sqr
	.p2align	4, 0x90
_fp_n_sqr:                              ## @fp_n_sqr
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	movl	%edx, %ebx
	movq	%rdi, %r14
	callq	_fp_sqr
	cmpl	$2, %ebx
	jl	LBB5_3
## %bb.1:
	decl	%ebx
	.p2align	4, 0x90
LBB5_2:                                 ## =>This Inner Loop Header: Depth=1
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	decl	%ebx
	jne	LBB5_2
LBB5_3:
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_from_mont                   ## -- Begin function fp_from_mont
	.p2align	4, 0x90
_fp_from_mont:                          ## @fp_from_mont
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%rdi, -56(%rbp)                 ## 8-byte Spill
	movq	(%rsi), %r9
	movq	8(%rsi), %rdi
	movq	16(%rsi), %rbx
	movq	24(%rsi), %rcx
	movq	%rcx, -48(%rbp)                 ## 8-byte Spill
	leaq	(%r9,%r9,4), %r8
	shlq	$56, %r8
	movabsq	$360287970189639680, %rax       ## imm = 0x500000000000000
	movq	%r9, %rdx
	mulxq	%rax, %r14, %r10
	movq	%rdi, %rdx
	mulxq	%rax, %rsi, %r11
	addq	%rcx, %r8
	addq	%r10, %rsi
	movq	%rbx, %rdx
	mulxq	%rax, %r10, %r12
	adcq	%r11, %r10
	movq	%r8, %rdx
	mulxq	%rax, %rdx, %r15
	adcq	%r12, %rdx
	adcq	$0, %r15
	xorl	%r11d, %r11d
	movq	%r9, %r12
	negq	%r12
	movl	$0, %r13d
	sbbq	%r13, %r13
	sarq	$63, %r13
	movq	%r13, %rcx
	subq	%rdi, %rcx
	sbbq	$0, %r13
	sarq	$63, %r13
	movq	%r13, %rax
	subq	%rbx, %rax
	sbbq	$0, %r13
	sarq	$63, %r13
	subq	%r8, %r14
	movl	$0, %r8d
	sbbq	%r8, %r8
	addq	%r13, %r14
	adcq	%r13, %r8
	sarq	$63, %r8
	addq	%r8, %rsi
	adcq	$0, %r8
	sarq	$63, %r8
	addq	%r8, %r10
	adcq	$0, %r8
	sarq	$63, %r8
	addq	%r8, %rdx
	adcq	$0, %r8
	sarq	$63, %r8
	addq	%r9, %r12
	adcq	%rdi, %rcx
	adcq	%rbx, %rax
	adcq	-48(%rbp), %r14                 ## 8-byte Folded Reload
	adcq	$0, %rsi
	adcq	$0, %r10
	adcq	$0, %rdx
	adcq	%r15, %r8
	movabsq	$-360287970189639680, %rax      ## imm = 0xFB00000000000000
	xorq	%r8, %rax
	movq	%rsi, %rcx
	andq	%r10, %rcx
	andq	%rdx, %rcx
	andq	%rax, %rcx
	addq	$1, %rcx
	adcq	$-1, %r11
	andq	%r11, %rsi
	movq	-56(%rbp), %rax                 ## 8-byte Reload
	movq	%rsi, (%rax)
	andq	%r11, %r10
	movq	%r10, 8(%rax)
	andq	%r11, %rdx
	movq	%rdx, 16(%rax)
	andq	%r8, %r11
	movq	%r11, 24(%rax)
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_to_mont                     ## -- Begin function fp_to_mont
	.p2align	4, 0x90
_fp_to_mont:                            ## @fp_to_mont
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	leaq	_R2(%rip), %rdx
	popq	%rbp
	jmp	_fp_mul                         ## TAILCALL
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_set                         ## -- Begin function fp_set
	.p2align	4, 0x90
_fp_set:                                ## @fp_set
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	(%rsi), %rax
	movq	%rax, (%rdi)
	movq	8(%rsi), %rax
	movq	%rax, 8(%rdi)
	movq	16(%rsi), %rax
	movq	%rax, 16(%rdi)
	movq	24(%rsi), %rax
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_set_zero                    ## -- Begin function fp_set_zero
	.p2align	4, 0x90
_fp_set_zero:                           ## @fp_set_zero
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	vxorps	%xmm0, %xmm0, %xmm0
	vmovups	%ymm0, (%rdi)
	popq	%rbp
	vzeroupper
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_set_one                     ## -- Begin function fp_set_one
	.p2align	4, 0x90
_fp_set_one:                            ## @fp_set_one
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	$51, (%rdi)
	vxorps	%xmm0, %xmm0, %xmm0
	vmovups	%xmm0, 8(%rdi)
	movabsq	$72057594037927936, %rax        ## imm = 0x100000000000000
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_is_zero                     ## -- Begin function fp_is_zero
	.p2align	4, 0x90
_fp_is_zero:                            ## @fp_is_zero
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	(%rdi), %rax
	movq	8(%rdi), %rcx
	movq	16(%rdi), %rdx
	movq	24(%rdi), %rsi
	movq	%rcx, %rdi
	orq	%rax, %rdi
	andq	%rax, %rcx
	andq	%rdx, %rcx
	orq	%rsi, %rdx
	orq	%rdi, %rdx
	notq	%rcx
	movabsq	$360287970189639679, %rdi       ## imm = 0x4FFFFFFFFFFFFFF
	xorq	%rsi, %rdi
	orq	%rcx, %rdi
	movq	%rdx, %rcx
	negq	%rcx
	orq	%rdx, %rcx
	movq	%rdi, %rax
	negq	%rax
	orq	%rdi, %rax
	andq	%rcx, %rax
	shrq	$63, %rax
	decl	%eax
                                        ## kill: def $eax killed $eax killed $rax
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_equal                       ## -- Begin function fp_equal
	.p2align	4, 0x90
_fp_equal:                              ## @fp_equal
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	xorl	%r8d, %r8d
	movq	(%rdi), %rcx
	subq	(%rsi), %rcx
	movq	8(%rdi), %rax
	movl	$0, %edx
	sbbq	%rdx, %rdx
	sarq	$63, %rdx
	subq	8(%rsi), %rax
	movl	$0, %r9d
	sbbq	%r9, %r9
	addq	%rdx, %rax
	adcq	%rdx, %r9
	sarq	$63, %r9
	movq	16(%rdi), %rdx
	subq	16(%rsi), %rdx
	movl	$0, %r10d
	sbbq	%r10, %r10
	addq	%r9, %rdx
	adcq	%r9, %r10
	sarq	$63, %r10
	movq	24(%rdi), %rdi
	subq	24(%rsi), %rdi
	movl	$0, %esi
	sbbq	%rsi, %rsi
	addq	%r10, %rdi
	adcq	%r10, %rsi
	sarq	$63, %rsi
	movl	%esi, %r9d
	andl	$2, %r9d
	subq	%r9, %rcx
	sbbq	%r8, %r8
	sarq	$63, %r8
	addq	%r8, %rax
	adcq	$0, %r8
	sarq	$63, %r8
	addq	%r8, %rdx
	adcq	$0, %r8
	shrq	$63, %r8
	movabsq	$-720575940379279360, %r9       ## imm = 0xF600000000000000
	andq	%rsi, %r9
	orq	%r8, %r9
	subq	%r9, %rdi
	movq	%rdi, %rsi
	shrq	$59, %rsi
	movq	%rsi, %r8
	negq	%r8
	movabsq	$-360287970189639680, %r9       ## imm = 0xFB00000000000000
	andq	%r8, %r9
	addq	%rcx, %rsi
	adcq	$0, %rax
	adcq	$0, %rdx
	adcq	%rdi, %r9
	movq	%rax, %rcx
	orq	%rsi, %rcx
	orq	%rdx, %rcx
	orq	%r9, %rcx
	andq	%rsi, %rax
	andq	%rdx, %rax
	notq	%rax
	movabsq	$360287970189639679, %rdx       ## imm = 0x4FFFFFFFFFFFFFF
	xorq	%r9, %rdx
	orq	%rax, %rdx
	movq	%rcx, %rsi
	negq	%rsi
	orq	%rcx, %rsi
	movq	%rdx, %rax
	negq	%rax
	orq	%rdx, %rax
	andq	%rsi, %rax
	shrq	$63, %rax
	decl	%eax
                                        ## kill: def $eax killed $eax killed $rax
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_sqrt                        ## -- Begin function fp_sqrt
	.p2align	4, 0x90
_fp_sqrt:                               ## @fp_sqrt
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	subq	$32, %rsp
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	movq	%rdi, %rbx
	vmovups	(%rdi), %ymm0
	vmovups	%ymm0, -48(%rbp)
	movq	%rdi, %rsi
	vzeroupper
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	leaq	-48(%rbp), %rdx
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_mul
	movl	$246, %r14d
	.p2align	4, 0x90
LBB13_1:                                ## =>This Inner Loop Header: Depth=1
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	decl	%r14d
	jne	LBB13_1
## %bb.2:
	addq	$32, %rsp
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_inv                         ## -- Begin function fp_inv
	.p2align	4, 0x90
_fp_inv:                                ## @fp_inv
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	subq	$32, %rsp
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	movq	%rdi, %rbx
	leaq	-48(%rbp), %r14
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_exp3div4
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%r14, %rsi
	movq	%rbx, %rdx
	callq	_fp_mul
	addq	$32, %rsp
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.p2align	4, 0x90                         ## -- Begin function fp_exp3div4
_fp_exp3div4:                           ## @fp_exp3div4
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r12
	pushq	%rbx
	subq	$160, %rsp
	.cfi_offset %rbx, -48
	.cfi_offset %r12, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%rsi, %r14
	movq	%rdi, %rbx
	leaq	-128(%rbp), %r15
	movq	%r15, %rdi
	callq	_fp_sqr
	leaq	-192(%rbp), %r12
	movq	%r12, %rdi
	movq	%r14, %rsi
	movq	%r15, %rdx
	callq	_fp_mul
	movq	%r15, %rdi
	movq	%r15, %rsi
	callq	_fp_sqr
	leaq	-160(%rbp), %r14
	movq	%r14, %rdi
	movq	%r12, %rsi
	movq	%r15, %rdx
	callq	_fp_mul
	leaq	-96(%rbp), %r15
	movq	%r15, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r15, %rdi
	movq	%r15, %rsi
	callq	_fp_sqr
	movq	%r15, %rdi
	movq	%r15, %rsi
	callq	_fp_sqr
	movq	%r15, %rdi
	movq	%r15, %rsi
	movq	%r14, %rdx
	callq	_fp_mul
	movq	%rbx, %rdi
	movq	%r15, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	movq	%r15, %rdx
	callq	_fp_mul
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	movq	%r14, %rdx
	callq	_fp_mul
	leaq	-64(%rbp), %r14
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	movq	%r14, %rdx
	callq	_fp_mul
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	movq	%r14, %rdx
	callq	_fp_mul
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	movq	%r14, %rdx
	callq	_fp_mul
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movl	$119, %r15d
	.p2align	4, 0x90
LBB15_1:                                ## =>This Inner Loop Header: Depth=1
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	decl	%r15d
	jne	LBB15_1
## %bb.2:
	leaq	-64(%rbp), %r14
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	movq	%r14, %rdx
	callq	_fp_mul
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	leaq	-96(%rbp), %rdx
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_mul
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%rbx, %rdi
	movq	%r14, %rsi
	movq	%rbx, %rdx
	callq	_fp_mul
	leaq	-128(%rbp), %rdx
	movq	%rbx, %rdi
	movq	%rbx, %rsi
	callq	_fp_mul
	addq	$160, %rsp
	popq	%rbx
	popq	%r12
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_fp_is_square                   ## -- Begin function fp_is_square
	.p2align	4, 0x90
_fp_is_square:                          ## @fp_is_square
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	subq	$32, %rsp
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	movq	%rdi, %rbx
	leaq	-48(%rbp), %r14
	movq	%r14, %rdi
	movq	%rbx, %rsi
	callq	_fp_exp3div4
	movq	%r14, %rdi
	movq	%r14, %rsi
	callq	_fp_sqr
	movq	%r14, %rdi
	movq	%r14, %rsi
	movq	%rbx, %rdx
	callq	_fp_mul
	movq	-48(%rbp), %rcx
	xorl	%esi, %esi
	addq	$-51, %rcx
	movl	$0, %edi
	adcq	$-1, %rdi
	sarq	$63, %rdi
	movq	-40(%rbp), %rax
	addq	%rdi, %rax
	adcq	$0, %rdi
	sarq	$63, %rdi
	movq	-32(%rbp), %rdx
	addq	%rdi, %rdx
	adcq	$0, %rdi
	sarq	$63, %rdi
	movq	-24(%rbp), %r9
	addq	%rdi, %r9
	adcq	$0, %rdi
	movabsq	$-72057594037927936, %r8        ## imm = 0xFF00000000000000
	addq	%r9, %r8
	adcq	$-1, %rdi
	sarq	$63, %rdi
	movl	%edi, %r9d
	andl	$2, %r9d
	subq	%r9, %rcx
	sbbq	%rsi, %rsi
	sarq	$63, %rsi
	addq	%rsi, %rax
	adcq	$0, %rsi
	sarq	$63, %rsi
	addq	%rsi, %rdx
	adcq	$0, %rsi
	shrq	$63, %rsi
	movabsq	$-720575940379279360, %r9       ## imm = 0xF600000000000000
	andq	%rdi, %r9
	orq	%rsi, %r9
	subq	%r9, %r8
	movq	%r8, %rsi
	shrq	$59, %rsi
	movq	%rsi, %rdi
	negq	%rdi
	movabsq	$-360287970189639680, %r9       ## imm = 0xFB00000000000000
	andq	%rdi, %r9
	addq	%rcx, %rsi
	adcq	$0, %rax
	adcq	$0, %rdx
	adcq	%r8, %r9
	movq	%rax, %rcx
	orq	%rsi, %rcx
	orq	%rdx, %rcx
	orq	%r9, %rcx
	andq	%rsi, %rax
	andq	%rdx, %rax
	notq	%rax
	movabsq	$360287970189639679, %rdx       ## imm = 0x4FFFFFFFFFFFFFF
	xorq	%r9, %rdx
	orq	%rax, %rdx
	movq	%rcx, %rax
	negq	%rax
	orq	%rcx, %rax
	movq	%rdx, %rcx
	negq	%rcx
	orq	%rdx, %rcx
	testq	%rax, %rcx
	setns	%al
	addq	$32, %rsp
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.section	__TEXT,__const
	.p2align	3, 0x0                          ## @R2
_R2:
	.quad	3689348814741912944             ## 0x3333333333333d70
	.quad	3689348814741910323             ## 0x3333333333333333
	.quad	3689348814741910323             ## 0x3333333333333333
	.quad	230584300921369395              ## 0x333333333333333

.subsections_via_symbols
