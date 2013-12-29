unsigned *ByteCodePtr;
Value_t *ImmedPtr;

ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;

Value_t x;
unsigned A;
Value_t y;
unsigned B;
unsigned C;
unsigned D;

switch(opcode) {
  
 TailCall_cAbs:
  
  if (!ByteCodePtr) {
    goto Laa;
  }
  
 case cAbs:
   switch( ByteCodePtr[0] ) {
   case cNeg:
     goto Lab;
   case cImmed:
     x=ImmedPtr[0];
     goto Lac;
   default: 
     A=ByteCodePtr[0];
     if (IsNeverNegativeValueOpcode(A)) {
       return;
     }
   }
   goto Default0;

 TailCall_cAdd:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cAdd :
   switch( ByteCodePtr[0] ) {
   case cDup: 
     switch( ByteCodePtr[-1] ) {
     case cAdd :
       if (ByteCodePtr[-2]==cDup) {
	 goto Lad;
       }
       break;
     case cMul :
       if (ByteCodePtr[-2]==cAdd) {
	 if (ByteCodePtr[-3]==cDup) {
	   goto Lae;
	 } 
       }
       break;
     } 
     goto Default1;
   case cMul :
     if (ByteCodePtr[-1]==cImmed) { 
       x=ImmedPtr[0];
       if (ByteCodePtr[-2]==cDup) {
	 goto Laf;
       } 
     } 
     goto Default1;
   case cNeg:
     goto Lag;
   case cImmed:
     x=ImmedPtr[0];
     switch( ByteCodePtr[-1] ) {
     case cAdd :
       if (ByteCodePtr[-2]==cImmed) { 
	 y=ImmedPtr[-1];
	 goto Lah;
       }
       break;
     case cRSub:
       switch( ByteCodePtr[-2] ) {
       case cAdd :
	 if (ByteCodePtr[-3]==cImmed) { 
	   y=ImmedPtr[-1];
	   goto Lai;
	 }
	 break;
       case cNeg:
	 goto Laj;
       case cImmed:
	 y=ImmedPtr[-1];
	 goto Lak;
       }
       break;
     case cSub :
       B=ByteCodePtr[-2];
       if (IsVarOpcode( B)) { 
	 switch( ByteCodePtr[-3] ) {
	 case cAdd :
	   if (ByteCodePtr[- 4 ]==cImmed) { 
	     y=ImmedPtr[-1];
	     goto Lal;
	   }
	   break;
	 case cNeg:
	   goto Lam;
	 case cImmed:
	   y=ImmedPtr[-1];
	   goto Lan;
	 } 
       }
       break;
     case cImmed:
       y=ImmedPtr[-1];
       goto Lao;
     } if (x==Value_t()) {
       goto Lap;
     }
     break;
   default: 
   Default1:
     A=ByteCodePtr[0];
     if (IsVarOpcode( A)) { 
       if (ByteCodePtr[-1]==cRSub ) {
	 if (ByteCodePtr[-2]==cImmed) { 
	   x=ImmedPtr[0];
	   incStackPtr();
	   --mStackPtr;
	   goto Laq;
	 };
	 incStackPtr();
	 --mStackPtr;
	 goto Lba;
       } 
     }
   }
   goto Default0;

 TailCall_cAnd:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cAnd :
   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lbb;
   case cImmed:
     x=ImmedPtr[0];
     if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lbc;
     }
     break;
   } 
   goto Default0;

 TailCall_cDiv:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cDiv :
   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lbd;
   case cImmed:
     x=ImmedPtr[0];
     if (x!=Value_t( )) {
       switch( ByteCodePtr[-1] ) {
       case cNeg:
	 goto Lbe;
       case cImmed:
	 y=ImmedPtr[-1];
	 goto Lbf;
       } 
     } if (x==Value_t( 1)) {
       goto Lap;
     }
     break;
   } 
   goto Default0;

 TailCall_cEqual:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cEqual:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 1)) { 
       A=ByteCodePtr[-1];
       if (IsLogicalOpcode(A)) {
	 goto Lap;
       } 
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lbg;
     } if (x==Value_t( 0 )) {
       switch( ByteCodePtr[-1] ) {
       case cAbs:
	 goto Lbh;
       case cSqr:
	 goto Lbh;
       };
       goto Lbi;
     }
   }
   goto Default0;

 TailCall_cGreater:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cGreater:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 0 )) {
       switch( ByteCodePtr[-1] ) {
       case cAbs:
	 goto Lbj;
       default: 
	 A=ByteCodePtr[-1];
	 if (IsNeverNegativeValueOpcode(A)) {
	   goto Lbk;
	 } 
       } 
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       ;
       goto Lbl;
     }
   }
   goto Default0;

 TailCall_cGreaterOrEq:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cGreaterOrEq:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 1 )) {
       switch( ByteCodePtr[-1] ) {
       case cAbs:
	 goto Lbj;
       default: 
	 A=ByteCodePtr[-1];
	 if (IsNeverNegativeValueOpcode(A)) {
	   goto Lbk;
	 } 
       } 
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lbm;
     }
   }
   goto Default0;

 TailCall_cInv:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cInv:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x!=Value_t()) {
       goto Lbn;
     }
   }
   goto Default0;
 TailCall_cLess:
   if (!ByteCodePtr) {
     goto Laa;
   }
 case cLess:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 0)) {A=ByteCodePtr[-1];
       if (IsNeverNegativeValueOpcode(A)) {
	 goto Lbo;
       } 
     } if (x==Value_t( 1 )) {
       switch( ByteCodePtr[-1] ) {
       case cAbs:
	 goto Lbp;
       default: 
	 A=ByteCodePtr[-1];
	 if (IsNeverNegativeValueOpcode(A)) {
	   goto Lbi;
	 } 
       } 
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lbq;
     }
   }
   goto Default0;

 TailCall_cLessOrEq:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cLessOrEq:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 0 )) {
       switch( ByteCodePtr[-1] ) {
       case cAbs:
	 goto Lbp;
       default: 
	 A=ByteCodePtr[-1];
	 if (IsNeverNegativeValueOpcode(A)) {
	   goto Lbi;
	 } 
       } 
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lca;
     }
   }
   goto Default0;

 TailCall_cMax:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cMax :
   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lcb;
   case cImmed:
     x=ImmedPtr[0];
     if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lcc;
     } 
     break;
   default: 
     A=ByteCodePtr[0];
     if (IsVarOpcode( A )) {
       switch( ByteCodePtr[-1] ) {
       case cDup:
	 B=ByteCodePtr[-2];
	 if (B==A) {
	   goto Lcb;
	 }
	 break;
       case cMax:
	 B=ByteCodePtr[-2];
	 if (B==A) {
	   goto Lcb;
	 }
	 break;
       } 
     }
   }
   goto Default0;

 TailCall_cMin:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cMin :
   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lcb;
   case cImmed:
     x=ImmedPtr[0];
     if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lcd;
     } 
     break;
   default: 
     A=ByteCodePtr[0];
     if (IsVarOpcode( A )) {
       switch( ByteCodePtr[-1] ) {
       case cDup:
	 B=ByteCodePtr[-2];
	 if (B==A) {
	   goto Lcb;
	 }
	 break;
       case cMin:
	 B=ByteCodePtr[-2];
	 if (B==A) {
	   goto Lcb;
	 }
	 break;
       } 
     }
   }
   goto Default0;

 TailCall_cMod:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cMod:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x!=Value_t()) {
       if (ByteCodePtr[-1]==cImmed) { 
	 y=ImmedPtr[-1];
	 goto Lce;
       } 
     }
   }
   goto Default0;

 TailCall_cMul:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cMul :
   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lcf;
   case cNeg: 
     switch( ByteCodePtr[-1] ) {
     case cDup:
       goto Lcg;
     default: 
       A=ByteCodePtr[-1];
       if (IsVarOpcode( A)) {
	 if (ByteCodePtr[-2]==cMul) { B=ByteCodePtr[-3];
	   if (B==A) {
	     goto Lch;
	   } 
	 } 
       } 
     } 
     goto Default2;
   case cPow :
     if (ByteCodePtr[-1]==cImmed) { 
       x=ImmedPtr[0];
       if (ByteCodePtr[-2]==cDup) {
	 goto Lci;
       } 
     } 
     goto Default2;
   case cImmed:
     x=ImmedPtr[0];
     if (x==Value_t( )) {
       switch( ByteCodePtr[-1] ) {
       case cMul:
	 A=ByteCodePtr[-2];
	 if (IsVarOpcode( A)) {
	   goto Lcj;
	 } 
	 goto Default3;
       default: Default3:
	 A=ByteCodePtr[-1];
	 if (IsBinaryOpcode(A)&&!HasInvalidRangesOpcode<
	     IsComplexType<Value_t>::result>( A)) { 
	   switch( ByteCodePtr[-2] ) {
	   case cImmed:
	     y=ImmedPtr[-1];
	     goto Lck;
	   default: B=ByteCodePtr[-2];
	     if (IsBinaryOpcode(B)&&!HasInvalidRangesOpcode<
		 IsComplexType<Value_t>::result>( B)) {
	       switch( ByteCodePtr[-3] ) {
	       case cImmed:
		 y=ImmedPtr[-1];
		 goto Lcl;
	       default: C=ByteCodePtr[-3];
		 if (IsVarOpcode( C)) {
		   goto Lcm;
		 } if (IsUnaryOpcode( C)&&!HasInvalidRangesOpcode<
		       IsComplexType<Value_t>::result>( C)) {
		   goto Lcn;
		 } 
	       } 
	     } if (IsVarOpcode( B)) {
	       goto Lcj;
	     } if (IsUnaryOpcode( B)&&!HasInvalidRangesOpcode<
		   IsComplexType<Value_t>::result>( B)) {
	       goto Lco;
	     } 
	   } 
	 } if (IsVarOpcode( A)) {
	   goto Lcp;
	 } if (IsUnaryOpcode( A)&&!HasInvalidRangesOpcode<
	       IsComplexType<Value_t>::result>( A)) {
	   goto Lcq;
	 } 
       } 
     } switch( ByteCodePtr[-1] ) {
     case cAdd:
       switch( ByteCodePtr[-2] ) {
       case cDup:
	 goto Lda;
       case cMul :
	 if (ByteCodePtr[-3]==cImmed) { 
	   y=ImmedPtr[-1];
	   A=ByteCodePtr[- 4];
	   if (IsVarOpcode( A)) {
	     goto Ldb;
	   } 
	 }
	 break;
       case cImmed:
	 y=ImmedPtr[-1];
	 goto Ldc;
       }
       break;
     case cMul:
       switch( ByteCodePtr[-2] ) {
       case cAdd :
	 if (ByteCodePtr[-3]==cDup) {
	   goto Ldd;
	 }
	 break;
       case cImmed:
	 y=ImmedPtr[-1];
	 if (y*x==Value_t( 1)) {
	   goto Lde;
	 };
	 goto Ldf;
       }
       break;
     case cNeg:
       goto Ldg;
     case cSub :
       if (ByteCodePtr[-2]==cMul) { 
	 if (ByteCodePtr[-3]==cImmed) { 
	   y=ImmedPtr[-1];
	   A=ByteCodePtr[- 4];
	   if (IsVarOpcode( A)) {
	     goto Ldh;
	   } 
	 } 
       }
       break;
     case cImmed:
       y=ImmedPtr[-1];
       goto Ldi;
     } if (x==Value_t( 1)) {
       goto Lap;
     } if (x==Value_t(-1)) {
       goto Ldj;
     } if (x==Value_t( 2)) {
       goto Ldk;
     }
     break;
   default: Default2:
     A=ByteCodePtr[0];
     if (IsVarOpcode( A )) {
       switch( ByteCodePtr[-1] ) {
       case cMul:
	 switch( ByteCodePtr[-2] ) {
	 case cNeg: B=ByteCodePtr[-3];
	   if (B==A) {
	     goto Lch;
	   } 
	   goto Default4;
	 default: Default4:
	   B=ByteCodePtr[-2];
	   if (B==A) {
	     goto Ldl;
	   } 
	 } 
	 goto Default5;
       case cNeg: B=ByteCodePtr[-2];
	 if (B==A) {
	   goto Lcg;
	 } 
	 goto Default5;
       default: Default5:
	 B=ByteCodePtr[-1];
	 if (B==A) {
	   goto Lcf;
	 } 
       } 
     } if (IsUnaryOpcode( A)) { 
       B=ByteCodePtr[-1];
       if (IsVarOpcode( B )&&mData->mByteCode .size()> 1 ) {
	 if (ByteCodePtr[-2]==cMul) { 
	   C=ByteCodePtr[-3];
	   if (C==A) {D=ByteCodePtr[- 4];
	     if (D==B) {
	       goto Ldm;
	     } 
	   } 
	 } 
       } 
     }
   }
   goto Default0;

 TailCall_cNEqual:

   if (!ByteCodePtr) {
     goto Laa;
   }

 case cNEqual:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 1)) { 
       A=ByteCodePtr[-1];
       if (IsLogicalOpcode(A)) {
	 goto Lbi;
       } 
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Ldn;
     } if (x==Value_t( 0 )) {
       switch( ByteCodePtr[-1] ) {
       case cAbs:
	 goto Ldo;
       case cSqr:
	 goto Ldo;
       };
       goto Lbk;
     }
   }
   goto Default0;

 TailCall_cNeg:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cNeg: 
   switch( ByteCodePtr[0] ) {
   case cMul :
     if (ByteCodePtr[-1]==cImmed) { 
       x=ImmedPtr[0];
       goto Ldp;
     }
     break;
   case cNeg:
     goto Lcb;
   case cImmed:
     x=ImmedPtr[0];
     goto Ldq;
   } 
   goto Default0;

 TailCall_cNot:

   if (!ByteCodePtr) {
     goto Laa;
   }

 case cNot:
   switch( ByteCodePtr[0] ) {
   case cAbs:
     goto Lea;
   case cAbsNot:
     A=ByteCodePtr[-1];
     if (IsLogicalOpcode(A)) {
       goto Lcb;
     } if (A!=cImmed) {
       goto Leb;
     } 
     goto Default6;
   case cAbsNotNot:
     goto Lec;
   case cAdd :
     if (ByteCodePtr[-1]==cImmed) { 
       x=ImmedPtr[0];
       goto Led;
     } 
     goto Default6;
   case cEqual:
     goto Lee;
   case cGreater:
     goto Lef;
   case cGreaterOrEq:
     goto Leg;
   case cLess:
     goto Leh;
   case cLessOrEq:
     goto Lei;
   case cNEqual:
     goto Lej;
   case cNeg:
     goto Lea;
   case cNot:
     goto Lbb;
   case cNotNot:
     goto Lea;
   case cImmed:
     x=ImmedPtr[0];
     goto Lek;
   default: 
   Default6:
     A=ByteCodePtr[0];
     if (IsNeverNegativeValueOpcode(A)) {
       goto Lel;
     }
   }
   goto Default0;

 TailCall_cNotNot:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cNotNot:
   switch( ByteCodePtr[0] ) {
   case cAdd :
     if (ByteCodePtr[-1]==cImmed) { 
       x=ImmedPtr[0];
       goto Lem;
     }
     break;
   case cNot:
     return;
   } 
   goto Default0;

 TailCall_cOr:

   if (!ByteCodePtr) {
     goto Laa;
   }

 case cOr:

   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lbb;
   case cImmed:
     x=ImmedPtr[0];
     if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Len;
     }
     break;
   } 
   goto Default0;

 TailCall_cRDiv:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cRDiv:
   if (ByteCodePtr[0 ]==cImmed) { 
     x=ImmedPtr[0];
     if (x==Value_t( 1)) {
       goto Leo;
     }
   }
   goto Default0;
 TailCall_cRSub:
   if (!ByteCodePtr) {
     goto Laa;
   }
 case cRSub :
   if (ByteCodePtr[0 ]==cDup) {
     goto Lep;
   } 
   goto Default0;

 TailCall_cSqr:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cSqr:
   switch( ByteCodePtr[0] ) {
   case cAbs:
     goto Leq;
   case cNeg:
     goto Leq;
   } 
   goto Default0;

 TailCall_cSub:

   if (!ByteCodePtr) {
     goto Laa;
   }
 case cSub :
   switch( ByteCodePtr[0] ) {
   case cDup:
     goto Lep;
   case cNeg:
     goto Lfa;
   case cImmed:
     x=ImmedPtr[0];
     if (x==Value_t()) {
       goto Lap;
     } if (ByteCodePtr[-1]==cImmed) { 
       y=ImmedPtr[-1];
       goto Lfb;
     };
     goto Lfc;
   default: 
     A=ByteCodePtr[0];
     if (IsVarOpcode( A)) { 
       if (ByteCodePtr[-1]==cRSub ) {
	 if (ByteCodePtr[-2]==cImmed) { 
	   x=ImmedPtr[0];
	   goto Lfd;
	 };
	 incStackPtr();
	 --mStackPtr;
	 goto Lfe;
       } 
     }
   }
   goto Default0;

 default: 

 Default0:

   A=opcode;
   if (IsComparisonOpcode(A)) { 
     if (ByteCodePtr[0 ]==cImmed) { 
       x=ImmedPtr[0];
       switch( ByteCodePtr[-1] ) {
       case cAdd :
	 if (ByteCodePtr[-2]==cImmed) { 
	   y=ImmedPtr[-1];
	   goto Lff;
	 }
	 break;
       case cNeg:
	 goto Lfg;
       } 
     } 
   } 
   if (IsVarOpcode( A )&&mData->mByteCode .size()> 0) {
     B=ByteCodePtr[0 ];
     if (B==A) {
       goto Lfh;
     } 
   } 
   if (IsUnaryOpcode( A)) { B=ByteCodePtr[0];
     if (IsVarOpcode( B )&&mData->mByteCode .size()> 1) {
       C=ByteCodePtr[-1];
       if (C==A) {D=ByteCodePtr[-2];
	 if (D==B) {
	   goto Lfi;
	 } 
       } 
     } 
   }
 } 
goto Laa;

Laa:
mData->mByteCode.push_back( opcode);
return;

Lab:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto TailCall_cAbs;

Lac:
ImmedPtr[0]=fp_abs(x);
return;

Lad:
mData->mImmed.push_back( Value_t( 4));
ByteCodePtr[-2]=cImmed;
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
ByteCodePtr-=2;

Lfj:
opcode=cMul;

Lfk:
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;

Lfl:
goto TailCall_cMul;

Lae:
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cMul);
mData->mImmed.push_back( Value_t( 4));
mData->mByteCode.push_back( cImmed );
goto Lfj;

Laf:
ImmedPtr[0]=x+Value_t( 1);
ByteCodePtr[-2]=cImmed;
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
ByteCodePtr-=2;

Lbo:
opcode=cMul;
goto Lfl;

Lag:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cSub;

Lfm:
goto TailCall_cSub;
Lah:ImmedPtr[-1]=y+x;
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();

Lfn:
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;

Lfo:
goto TailCall_cAdd;

Lai:
ImmedPtr[-1]=y+x;
mData->mImmed.pop_back();
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();

Lfp:
AddFunctionOpcode( cAdd);

Lfq:
opcode=cRSub;
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;
goto TailCall_cRSub;

Laj:
ImmedPtr[0]=-x;
ByteCodePtr[-2]=cImmed;
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
goto Lfp;

Lak:
ImmedPtr[-1]=y+x;
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
goto Lfq;

Lal:
ImmedPtr[-1]=y+x;
mData->mImmed.pop_back();
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();

Lga:
AddFunctionOpcode( cAdd);

Lgb:
AddFunctionOpcode( B);

Lgc:
opcode=cSub;
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;
goto Lfm;

Lam:
ImmedPtr[0]=-x;
ByteCodePtr[-3]=cImmed;
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
goto Lga;

Lan:
ImmedPtr[-1]=y+x;
mData->mImmed.pop_back();
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
goto Lgb;

Lao:
ImmedPtr[-1]=y+x;

Lap:
mData->mImmed.pop_back();

Lcb:
mData->mByteCode.pop_back();
return;

Laq:
mData->mImmed.pop_back();
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( A);
mData->mImmed.push_back(x );
mData->mByteCode.push_back( cImmed );
goto Lfp;

Lba:
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( A );
goto Lfp;

Lbb:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cNotNot;

Lgd:
goto TailCall_cNotNot;

Lbc:
ImmedPtr[-1]=fp_and(x ,y );
goto Lap;

Lbd:
mData->mImmed.push_back( Value_t());
ByteCodePtr[0]=cImmed;
AddFunctionOpcode( cMul);
mData->mImmed.push_back( Value_t( 1));

Lge:
mData->mByteCode.push_back( cImmed);

Lgf:
opcode=cAdd;
goto Lfn;

Lbe:
ImmedPtr[0]=-x;
ByteCodePtr[-1]=cImmed;
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto TailCall_cDiv;

Lbf:
ImmedPtr[-1]=y/x;
goto Lap;

Lbg:
ImmedPtr[-1]=fp_equal (y ,x );
goto Lap;

Lbh:
ByteCodePtr[-1]=cImmed;
mData->mByteCode.pop_back();
ByteCodePtr-=1;

Lgg:
goto TailCall_cEqual;

Lbi:
mData->mImmed.pop_back();
mData->mByteCode.pop_back();

Lgh:
opcode=cNot;
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;

Lgi:
goto TailCall_cNot;

Lbj:
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();

Lgj:
opcode=cNotNot;
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;
goto Lgd;

Lbk:
mData->mImmed.pop_back();
mData->mByteCode.pop_back();
goto Lgj;

Lbl:
ImmedPtr[-1]=fp_less(x ,y );
goto Lap;

Lbm:
ImmedPtr[-1]=fp_lessOrEq(x ,y );
goto Lap;

Lbn:
ImmedPtr[0]=Value_t( 1)/x;
return;

Lbp:
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
goto Lgh;

Lbq:
ImmedPtr[-1]=fp_less (y ,x );
goto Lap;

Lca:
ImmedPtr[-1]=fp_lessOrEq (y ,x );
goto Lap;

Lcc:
ImmedPtr[-1]=fp_max(x ,y );
goto Lap;

Lcd:
ImmedPtr[-1]=fp_min(x ,y );
goto Lap;

Lce:
ImmedPtr[-1]=fp_mod (y ,x );
goto Lap;

Lcf:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cSqr;

Lgk:
goto TailCall_cSqr;

Lcg:
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cSqr);

Lgl:
opcode=cNeg;
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;
goto TailCall_cNeg;

Lch:
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cSqr);
AddFunctionOpcode( cMul );
goto Lgl;

Lci:
ImmedPtr[0]=x+Value_t( 1);
ByteCodePtr[-2]=cImmed;
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cPow);
return;

Lcj:
ByteCodePtr[-2]=cImmed;
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
ByteCodePtr-=2;
goto Lfl;

Lck:
ImmedPtr[-1]=x;

Lgm:
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
goto Lfk;

Lcl:
for (unsigned tmp=2;tmp-->0;) mData->mImmed.pop_back();

Lgn:
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();

Lgo:
AddFunctionOpcode( A );
mData->mImmed.push_back( x);

Lgp:
mData->mByteCode.push_back( cImmed );
goto Lfk;

Lcm:
mData->mImmed.pop_back();
goto Lgn;

Lcn:
mData->mImmed.pop_back();
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( B );
goto Lgo;

Lco:
mData->mImmed.pop_back();
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
goto Lgo;

Lcp:
ByteCodePtr[-1]=cImmed;
goto Lcb;

Lcq:
ByteCodePtr[-1]=cImmed;
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto Lfl;

Lda:
ImmedPtr[0]=x+x;
goto Lcj;

Ldb:
ImmedPtr[-1]=x;
ByteCodePtr[- 4]=cImmed;
mData->mImmed.pop_back();
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cMul);
AddFunctionOpcode( A );
mData->mImmed.push_back( y*x );
mData->mByteCode.push_back( cImmed);
AddFunctionOpcode( cMul );
goto Lgf;

Ldc:
ImmedPtr[-1]=x;
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cMul);
mData->mImmed.push_back( y*x );
goto Lge;

Ldd:
mData->mImmed.pop_back();
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cMul);
mData->mImmed.push_back( x+x );
goto Lgp;

Lde:
for (unsigned tmp=2;tmp-->0;) mData->mImmed.pop_back();
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
return;

Ldf:
ImmedPtr[-1]=y*x;
goto Lgm;

Ldg:
ImmedPtr[0]=-x;
goto Lcq;

Ldh:
ImmedPtr[-1]=x;
ByteCodePtr[- 4]=cImmed;
mData->mImmed.pop_back();
for (unsigned tmp=4;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( cMul);
AddFunctionOpcode( A );
mData->mImmed.push_back( y*x );
mData->mByteCode.push_back( cImmed);
AddFunctionOpcode( cMul );
goto Lgc;

Ldi:
ImmedPtr[-1]=y*x;
goto Lap;

Ldj:
mData->mImmed.pop_back();
mData->mByteCode.pop_back();
goto Lgl;

Ldk:
ByteCodePtr[0]=cDup;
ImmedPtr-=1;
mData->mImmed.pop_back();

Lgq:
opcode=cAdd;
goto Lfo;

Ldl:
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();

Lha:
AddFunctionOpcode( cSqr );
goto Lfk;

Ldm:
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
goto Lha;

Ldn:
ImmedPtr[-1]=fp_nequal (y ,x );
goto Lap;

Ldo:
ByteCodePtr[-1]=cImmed;
mData->mByteCode.pop_back();
ByteCodePtr-=1;

Lhb:
goto TailCall_cNEqual;

Ldp:
ImmedPtr[0]=-x;
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto Lbo;

Ldq:
ImmedPtr[0]=-x;
return;

Lea:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto Lgi;

Leb:
mData->mByteCode.pop_back();
AddFunctionOpcode( cAbsNotNot);
return;

Lec:
mData->mByteCode.pop_back();
Lel:AddFunctionOpcode( cAbsNot);
return;

Led:
ImmedPtr[0]=-x;

Lej:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cEqual;
goto Lgg;

Lee:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cNEqual;
goto Lhb;

Lef:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cLessOrEq;
goto TailCall_cLessOrEq;

Leg:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cLess;
goto TailCall_cLess;

Leh:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cGreaterOrEq;
goto TailCall_cGreaterOrEq;

Lei:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
opcode=cGreater;
goto TailCall_cGreater;

Lek:
ImmedPtr[0]=fp_not (x);
return;

Lem:
ImmedPtr[0]=-x;
goto Lee;

Len:
ImmedPtr[-1]=fp_or(x ,y );
goto Lap;

Leo:
mData->mImmed.pop_back();
mData->mByteCode.pop_back();
opcode=cInv;
ByteCodePtr=!mData->mByteCode.empty() ? 
  &mData->mByteCode[0]+mData->mByteCode .size()-1:0;
ImmedPtr=!mData->mImmed.empty() ? 
  &mData->mImmed[0]+mData->mImmed .size()-1:0;
goto TailCall_cInv;

Lep:
mData->mImmed.push_back( Value_t());
ByteCodePtr[0]=cImmed;
goto Lfj;

Leq:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto Lgk;

Lfa:
mData->mByteCode.pop_back();
ByteCodePtr-=1;
goto Lgq;

Lfb:
ImmedPtr[-1]=y-x;
goto Lap;

Lfc:
ImmedPtr[0]=-x;
goto Lgq;

Lfd:
mData->mImmed.pop_back();
for (unsigned tmp=3;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( A );
AddFunctionOpcode( cAdd );
mData->mImmed.push_back( x );
mData->mByteCode.push_back( cImmed );
goto Lfq;

Lfe:
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( A );
AddFunctionOpcode( cSub );
goto Lfq;

Lff:
ImmedPtr[-1]=x-y;
mData->mImmed.pop_back();
for (unsigned tmp=2;tmp-->0;) mData->mByteCode.pop_back();
AddFunctionOpcode( A);
return;

Lfg:
ImmedPtr[0]=-x;
ByteCodePtr[-1]=cImmed;
mData->mByteCode.pop_back();
AddFunctionOpcode( OppositeComparisonOpcode(A));
return;

Lfh:
mData->mByteCode.push_back( cDup);
return;

Lfi:
ByteCodePtr[0]=cDup;
return;
