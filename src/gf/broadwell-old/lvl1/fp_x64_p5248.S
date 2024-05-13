//*******************************************************************************************
// Field arithmetic in x64 assembly for 5·2²⁴⁸ - 1 on Linux
//*******************************************************************************************  

.intel_syntax noprefix 

#define _ADX_
  
// Format function and variable names for Mac OS X
#if defined(__APPLE__)
    #define fmt(f)    _##f
#else
    #define fmt(f)    f
#endif

// Registers that are used for parameter passing:
#define reg_p1  rdi
#define reg_p2  rsi
#define reg_p3  rdx

.p_lows: .quad 0xffffffffffffffff
.p_high: .quad 0x04ffffffffffffff 	

.text
.p2align 4,,15

//***********************************************************************
//  Field addition
//  Operation: c [reg_p1] = a [reg_p2] + b [reg_p3]
//*********************************************************************** 
.global fmt(fpadd5248_asm)
fmt(fpadd5248_asm):
  xor    rax, rax
  mov    r8, [reg_p2]
  mov    r9, [reg_p2+8]
  mov    r10, [reg_p2+16]
  mov    r11, [reg_p2+24]
  add    r8, [reg_p3] 
  adc    r9, [reg_p3+8] 
  adc    r10, [reg_p3+16] 
  adc    r11, [reg_p3+24] 

  mov    rcx, [rip + .p_lows]
  sub    r8, rcx
  sbb    r9, rcx
  sbb    r10, rcx
  mov    rsi, [rip + .p_high]
  sbb    r11, rsi
  sbb    rax, 0
  
  and    rcx, rax
  and    rsi, rax
  
  add    r8, rcx
  adc    r9, rcx  
  adc    r10, rcx  
  adc    r11, rsi 
  mov    [reg_p1], r8
  mov    [reg_p1+8], r9 
  mov    [reg_p1+16], r10 
  mov    [reg_p1+24], r11
  
  ret


//***********************************************************************
//  Field subtraction
//  Operation: c [reg_p1] = a [reg_p2] - b [reg_p3]
//*********************************************************************** 
.global fmt(fpsub5248_asm)
fmt(fpsub5248_asm):
  xor    rax, rax
  mov    r8, [reg_p2]
  mov    r9, [reg_p2+8]
  mov    r10, [reg_p2+16]
  mov    r11, [reg_p2+24]
  sub    r8, [reg_p3] 
  sbb    r9, [reg_p3+8] 
  sbb    r10, [reg_p3+16] 
  sbb    r11, [reg_p3+24] 
  sbb    rax, 0
  
  mov    rcx, [rip + .p_lows]
  mov    rsi, [rip + .p_high]
  and    rcx, rax
  and    rsi, rax
  add    r8, rcx  
  adc    r9, rcx
  adc    r10, rcx  
  adc    r11, rsi 
  mov    [reg_p1], r8
  mov    [reg_p1+8], r9 
  mov    [reg_p1+16], r10 
  mov    [reg_p1+24], r11 
  ret


///////////////////////////////////////////////////////////////// MACRO
.macro SUB5248_PX
  mov    r8, [reg_p1]
  mov    r9, [reg_p1+8]
  mov    r10, [reg_p1+16]
  mov    r11, [reg_p1+24]
  sub    r8, [reg_p2] 
  sbb    r9, [reg_p2+8] 
  sbb    r10, [reg_p2+16] 
  sbb    r11, [reg_p2+24] 

  mov    rax, [rip + .p_lows]
  mov    rsi, [rip + .p_high]
  add    r8, rax
  adc    r9, rax
  adc    r10, rax 
  adc    r11, rsi 
  mov    [reg_p3], r8
  mov    [reg_p3+8], r9 
  mov    [reg_p3+16], r10 
  mov    [reg_p3+24], r11
  .endm


//***********************************************************************
//  Multiprecision subtraction with correction with p5248
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2] + p5248
//*********************************************************************** 
.global fmt(mp_sub5248_p2_asm)
fmt(mp_sub5248_p2_asm):

  SUB5248_PX
  ret


///////////////////////////////////////////////////////////////// MACRO
// Schoolbook integer multiplication
// Inputs:  memory pointers M0 and M1
// Outputs: memory pointer C and regs T1, T3, rax
// Temps:   regs T0:T6
/////////////////////////////////////////////////////////////////

#ifdef _ADX_
.macro MUL192_SCHOOL M0, M1, C, T0, T1, T2, T3, T4, T5, T6
    mov    rdx, \M0
    mulx   \T0, \T1, \M1     // T0:T1 = A0*B0
    mov    \C, \T1           // C0_final
    mulx   \T1, \T2, 8\M1    // T1:T2 = A0*B1
    xor    rax, rax   
    adox   \T0, \T2        
    mulx   \T2, \T3, 16\M1   // T2:T3 = A0*B2
    adox   \T1, \T3
           
    mov    rdx, 8\M0
    mulx   \T3, \T4, \M1     // T3:T4 = A1*B0
    adox   \T2, rax 
    xor    rax, rax   
    mulx   \T5, \T6, 8\M1    // T5:T6 = A1*B1
    adox   \T4, \T0
    mov    8\C, \T4          // C1_final  
    adcx   \T3, \T6      
    mulx   \T6, \T0, 16\M1   // T6:T0 = A1*B2 
    adox   \T3, \T1  
    adcx   \T5, \T0     
    adcx   \T6, rax 
    adox   \T5, \T2	
    
    mov    rdx, 16\M0
    mulx   \T1, \T0, \M1     // T1:T0 = A2*B0
    adox   \T6, rax
    xor    rax, rax 
    mulx   \T4, \T2, 8\M1    // T4:T2 = A2*B1
    adox   \T0, \T3   
    mov    16\C, \T0         // C2_final 
    adcx   \T1, \T5    
    mulx   \T0, \T3, 16\M1   // T0:T3 = A2*B2
    adcx   \T4, \T6  
    adcx   \T0, rax
    adox   \T1, \T2
    adox   \T3, \T4
    adox   rax, \T0
.endm 
    
///////////////////////////////////////////////////////////////// MACRO
// Schoolbook integer multiplication
// Inputs:  memory pointers M0 and M1
// Outputs: memory pointer C
// Temps:   regs T0:T9
/////////////////////////////////////////////////////////////////

.macro MUL256_SCHOOL M0, M1, C, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9
    mov    rdx, \M0
    mulx   \T0, \T1, \M1     // T0:T1 = A0*B0
    mov    \C, \T1           // C0_final
    mulx   \T1, \T2, 8\M1    // T1:T2 = A0*B1
    xor    rax, rax   
    adox   \T0, \T2        
    mulx   \T2, \T3, 16\M1   // T2:T3 = A0*B2
    adox   \T1, \T3        
    mulx   \T3, \T4, 24\M1   // T3:T4 = A0*B3
    adox   \T2, \T4 
           
    mov    rdx, 8\M0
    mulx   \T5, \T4, \M1     // T5:T4 = A1*B0
    adox   \T3, rax 
    xor    rax, rax   
    mulx   \T6, \T7, 8\M1    // T6:T7 = A1*B1
    adox   \T4, \T0
    mov    8\C, \T4          // C1_final  
    adcx   \T5, \T7      
    mulx   \T7, \T8, 16\M1   // T7:T8 = A1*B2
    adcx   \T6, \T8  
    adox   \T5, \T1      
    mulx   \T8, \T9, 24\M1   // T8:T9 = A1*B3
    adcx   \T7, \T9        
    adcx   \T8, rax   
    adox   \T6, \T2
    
    mov    rdx, 16\M0
    mulx   \T1, \T0, \M1     // T1:T0 = A2*B0
    adox   \T7, \T3
    adox   \T8, rax
    xor    rax, rax 
    mulx   \T2, \T3, 8\M1    // T2:T3 = A2*B1
    adox   \T0, \T5   
    mov    16\C, \T0         // C2_final 
    adcx   \T1, \T3    
    mulx   \T3, \T4, 16\M1   // T3:T4 = A2*B2
    adcx   \T2, \T4 
    adox   \T1, \T6       
    mulx   \T4,\T9, 24\M1    // T3:T4 = A2*B3
    adcx   \T3, \T9        
    mov    rdx, 24\M0
    adcx   \T4, rax         

    adox   \T2, \T7
    adox   \T3, \T8
    adox   \T4, rax

    mulx   \T5, \T0, \M1     // T5:T0 = A3*B0
    xor    rax, rax 
    mulx   \T6, \T7, 8\M1    // T6:T7 = A3*B1
    adcx   \T5, \T7 
    adox   \T1, \T0       
    mulx   \T7, \T8, 16\M1   // T7:T8 = A3*B2
    adcx   \T6, \T8  
    adox   \T2, \T5      
    mulx   \T8, \T9, 24\M1   // T8:T9 = A3*B3
    adcx   \T7, \T9        
    adcx   \T8, rax         

    adox   \T3, \T6
    adox   \T4, \T7
    adox   \T8, rax
    mov    24\C, \T1         // C3_final
    mov    32\C, \T2         // C4_final
    mov    40\C, \T3         // C5_final
    mov    48\C, \T4         // C6_final
    mov    56\C, \T8         // C7_final
.endm 

#else

.macro MUL192_SCHOOL M0, M1, C, T0, T1, T2, T3, T4, T5, T6
    mov    rdx, \M0
    mulx   \T0, \T1, \M1     // T0:T1 = A0*B0
    mov    \C, \T1           // C0_final
    mulx   \T1, \T2, 8\M1    // T1:T2 = A0*B1
    add    \T0, \T2        
    mulx   \T2, \T3, 16\M1   // T2:T3 = A0*B2
    adc    \T1, \T3
           
    mov    rdx, 8\M0
    mulx   \T3, \T4, \M1     // T3:T4 = A1*B0
    adc    \T2, 0   
    mulx   \T5, \T6, 8\M1    // T5:T6 = A1*B1
    add    \T4, \T0
    mov    8\C, \T4          // C1_final
    adc    \T3, \T1  
    adc    \T5, \T2	    
    mulx   \T2, \T1, 16\M1   // T2:T1 = A1*B2
    adc    \T2, 0    

    add    \T3, \T6  
    adc    \T5, \T1     
    adc    \T2, 0
    
    mov    rdx, 16\M0
    mulx   \T1, \T0, \M1     // T1:T0 = A2*B0
    add    \T0, \T3   
    mov    16\C, \T0         // C2_final 
    mulx   \T4, \T6, 8\M1    // T4:T6 = A2*B1
    adc    \T1, \T5    
    adc    \T2, \T4 
    mulx   rax, \T3, 16\M1   // rax:T3 = A2*B2 
    adc    rax, 0
    add    \T1, \T6
    adc    \T3, \T2
    adc    rax, 0
.endm 

.macro MUL256_SCHOOL M0, M1, C, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9
    mov    rdx, \M0
    mulx   \T0, \T1, \M1     //          T0:T1 = A0*B0
    mov    \C, \T1           // C0_final
    mulx   \T1, \T2, 8\M1    //       T1:T2 = A0*B1
    add    \T0, \T2        
    mulx   \T2, \T3, 16\M1   //    T2:T3 = A0*B2
    adc    \T1, \T3         
    mulx   \T3, \T4, 24\M1   // T3:T4 = A0*B3
    adc    \T2, \T4        
    mov    rdx, 8\M0
    adc    \T3, 0         

    mulx   \T5, \T4, \M1     //          T5:T4 = A1*B0
    mulx   \T6, \T7, 8\M1    //       T6:T7 = A1*B1
    add    \T5, \T7        
    mulx   \T7, \T8, 16\M1   //    T7:T8 = A1*B2
    adc    \T6, \T8        
    mulx   \T8, \T9, 24\M1   // T8:T9 = A1*B3
    adc    \T7, \T9        
    adc    \T8, 0         

    add    \T4, \T0
    mov    8\C, \T4          // C1_final
    adc    \T5, \T1
    adc    \T6, \T2
    adc    \T7, \T3
    mov    rdx, 16\M0
    adc    \T8, 0

    mulx   \T1, \T0, \M1     //          T1:T0 = A2*B0
    mulx   \T2, \T3, 8\M1    //       T2:T3 = A2*B1
    add    \T1, \T3        
    mulx   \T3, \T4, 16\M1   //    T3:T4 = A2*B2
    adc    \T2, \T4        
    mulx   \T4, \T9, 24\M1   // T4:T9 = A2*B3
    adc    \T3, \T9        
    mov    rdx, 24\M0
    adc    \T4, 0          

    add    \T0, \T5
    mov    16\C, \T0         // C2_final
    adc    \T1, \T6
    adc    \T2, \T7
    adc    \T3, \T8
    adc    \T4, 0

    mulx   \T5, \T0, \M1     //          T5:T0 = A3*B0
    mulx   \T6, \T7, 8\M1    //       T6:T7 = A3*B1
    add    \T5, \T7        
    mulx   \T7, \T8, 16\M1   //    T7:T8 = A3*B2
    adc    \T6, \T8        
    mulx   \T8, \T9, 24\M1   // T8:T9 = A3*B3
    adc    \T7, \T9         
    adc    \T8, 0         

    add    \T1, \T0
    mov    24\C, \T1         // C3_final
    adc    \T2, \T5
    mov    32\C, \T2         // C4_final
    adc    \T3, \T6
    mov    40\C, \T3         // C5_final
    adc    \T4, \T7
    mov    48\C, \T4         // C6_final
    adc    \T8, 0
    mov    56\C, \T8         // C7_final
.endm
#endif


//*****************************************************************************
//  256-bit multiplication using schoolbook
//***************************************************************************** 
.global fmt(mul5248_asm)
fmt(mul5248_asm):    
    push   r12
    push   r13 
    push   r14 
    push   r15
    push   rbx 
    push   rbp
    mov rcx, reg_p3

    MUL256_SCHOOL  [reg_p2], [rcx], [reg_p1], r8, r9, r10, r11, r12, r13, r14, r15, rbx, rbp

    pop    rbp
    pop    rbx
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret
  
//**************************************************************************************
//  Montgomery reduction
//  Operation: c [reg_p1] = a [reg_p2]
//************************************************************************************** 
.global fmt(rdc5248_asm)
fmt(rdc5248_asm):
  push   r12
  push   r13 
  push   r14 
  push   r15
  
  xor rax, rax
  // r8  = a[0]
  // r11 = (a[0] << 58) + (a[0] << 56)
  // rax = carry
  mov r8, [reg_p2]
	mov r11, r8
  shl r11, 58
  mov r9, r8
  shl r9, 56
  add r11, r9
  adc rax, 0
  // r11 += a[3]
  // rax += carry
  add r11, [reg_p2+24]
  adc rax, 0

  // [r9:r10] = a[1:2]
  // [r12:r15] = [r8:r11]
  mov r9, [reg_p2+8]
  mov r10, [reg_p2+16]
  mov r12, r8
  mov r13, r9
  mov r14, r10
  mov r15, r11

  // [r8:r11] >>= 6
  // [r12:r15] >>= 8
  shrd r8, r9, 6
  shrd r12, r13, 8
  shrd r9, r10, 6
  shrd r13, r14, 8
  shrd r10, r11, 6
  shrd r14, r15, 8
  shr r11, 6
  shr r15, 8

  // [r12:r15] += rax + [r12:r15]
  or rax, rax
  adox r8, rax
  adox r9, r13
  adox r10, r14
  adox r11, r15
  xor rax, rax
  adcx r8, r12
  adcx r9, rax
  adcx r10, rax
  adcx r11, rax
  // [r8:r11] += a[4-7]
  xor rax, rax
  adox r8, [reg_p2+32]
  adox r9, [reg_p2+40]
  adox r10, [reg_p2+48]
  adox r11, [reg_p2+56]

  // remove p5248 if > p5248
  mov    rcx, [rip + .p_lows]
  sub    r8, rcx
  sbb    r9, rcx
  sbb    r10, rcx
  mov    rdx, [rip + .p_high]
  sbb    r11, rdx
  sbb    rax, 0
  
  and    rcx, rax
  and    rdx, rax
  
  add    r8, rcx
  adc    r9, rcx
  adc    r10, rcx  
  adc    r11, rdx 

  mov    [reg_p1], r8
  mov    [reg_p1+8], r9 
  mov    [reg_p1+16], r10 
  mov    [reg_p1+24], r11

  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret

//***********************************************************************
//  255-bit multiprecision addition
//  Operation: c [reg_p1] = a [reg_p2] + b [reg_p3]
//*********************************************************************** 
.global fmt(mp_add255_asm)
fmt(mp_add255_asm): 
  mov    r8, [reg_p2]
  mov    r9, [reg_p2+8]
  mov    r10, [reg_p2+16]
  mov    r11, [reg_p2+24]
  add    r8, [reg_p3] 
  adc    r9, [reg_p3+8] 
  adc    r10, [reg_p3+16] 
  adc    r11, [reg_p3+24] 
  mov    [reg_p1], r8
  mov    [reg_p1+8], r9
  mov    [reg_p1+16], r10
  mov    [reg_p1+24], r11
  ret


//***************************************************************************
//  2x255-bit multiprecision subtraction/addition
//  Operation: c [reg_p1] = a [reg_p2] - b [reg_p3]. If c < 0, add p5248*2^256
//*************************************************************************** 
.global fmt(mp_subadd255x2_asm)
fmt(mp_subadd255x2_asm):
  push   r12
  push   r13 
  xor    rax, rax

  mov    r8, [reg_p2]
  mov    r9, [reg_p2+8]
  mov    r10, [reg_p2+16]
  mov    r11, [reg_p2+24]
  sub    r8, [reg_p3] 
  sbb    r9, [reg_p3+8] 
  sbb    r10, [reg_p3+16] 
  sbb    r11, [reg_p3+24] 
  mov    [reg_p1], r8
  mov    [reg_p1+8], r9
  mov    [reg_p1+16], r10
  mov    [reg_p1+24], r11

  mov    r8, [reg_p2+32]
  mov    r9, [reg_p2+40]
  mov    r10, [reg_p2+48] 
  mov    r11, [reg_p2+56]
  sbb    r8, [reg_p3+32] 
  sbb    r9, [reg_p3+40] 
  sbb    r10, [reg_p3+48]
  sbb    r11, [reg_p3+56] 
  sbb    rax, 0
  
  // Add p5248 anded with the mask in rax 
  mov    r12, [rip + .p_lows]
  mov    r13, [rip + .p_high]
  and    r12, rax
  and    r13, rax
  add    r8, r12
  adc    r9, r12
  adc    r10, r12
  adc    r11, r13
  
  mov    [reg_p1+32], r8
  mov    [reg_p1+40], r9
  mov    [reg_p1+48], r10
  mov    [reg_p1+56], r11
  pop    r13
  pop    r12
  ret


//***********************************************************************
//  Double 2x255-bit multiprecision subtraction
//  Operation: c [reg_p1] = c [reg_p1] - a [reg_p2] - b [reg_p3]
//*********************************************************************** 
.global fmt(mp_dblsub255x2_asm)
fmt(mp_dblsub255x2_asm):
  push   r12
  push   r13
  push   r14
  push   r15
  
  mov    r8, [reg_p1]
  mov    r9, [reg_p1+8]
  mov    r10, [reg_p1+16]
  mov    r11, [reg_p1+24]
  mov    r12, [reg_p1+32]
  mov    r13, [reg_p1+40]
  mov    r14, [reg_p1+48]
  mov    r15, [reg_p1+48]
  sub    r8, [reg_p2]
  sbb    r9, [reg_p2+8] 
  sbb    r10, [reg_p2+16] 
  sbb    r11, [reg_p2+24] 
  sbb    r12, [reg_p2+32] 
  sbb    r13, [reg_p2+40] 
  sbb    r14, [reg_p2+48]
  sbb    r15, [reg_p2+56]
  sub    r8, [reg_p3]
  sbb    r9, [reg_p3+8] 
  sbb    r10, [reg_p3+16] 
  sbb    r11, [reg_p3+24] 
  sbb    r12, [reg_p3+32] 
  sbb    r13, [reg_p3+40] 
  sbb    r14, [reg_p3+48]
  sbb    r15, [reg_p3+56]
  mov    [reg_p1], r8
  mov    [reg_p1+8], r9
  mov    [reg_p1+16], r10
  mov    [reg_p1+24], r11
  mov    [reg_p1+32], r12
  mov    [reg_p1+40], r13
  mov    [reg_p1+48], r14
  mov    [reg_p1+56], r15

  pop    r15
  pop    r14
  pop    r13
  pop    r12
  ret
