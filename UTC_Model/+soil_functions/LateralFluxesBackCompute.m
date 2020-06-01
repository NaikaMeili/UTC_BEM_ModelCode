function[SQ2to1, SQ1to2, SQ3to1, SQ3to2, SQ2to3, SQ1to3]=LateralFluxesBackCompute(Q1,Q2,Q3,f1,f2,f3)

syms Q2to1 Q3to1 Q1to2 Q3to2 Q2to3 Q1to3

eqn1	=	Q2to1 + Q3to1 == Q1;
eqn2	=	Q1to2 + Q3to2 == Q2;
eqn3	=	Q2to3 + Q1to3 == Q3;
eqn4	=	Q2to1*f1- Q1to2*f2 == 0;
eqn5	=	Q3to1*f1 - Q1to3*f3 == 0;	
eqn6	=	Q3to2*f2 - Q2to3*f3 == 0;

[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [Q2to1 Q3to1 Q1to2 Q3to2 Q2to3 Q1to3]);

X = linsolve(A,B);

SQ2to1	=	double(X(1));
SQ3to1	=	double(X(2));
SQ1to2	=	double(X(3));
SQ3to2	=	double(X(4));
SQ2to3	=	double(X(5));
SQ1to3	=	double(X(6));


end


