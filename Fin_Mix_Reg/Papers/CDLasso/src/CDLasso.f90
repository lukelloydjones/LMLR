      MODULE COORDINATE_DESCENT
!
!     Determine double precision and set constants.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBLE = KIND(0.0D0)
      REAL(KIND=DBLE), PARAMETER :: ZERO  = 0.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: ONE   = 1.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: TWO   = 2.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: THREE = 3.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: FOUR  = 4.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: FIVE  = 5.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: SIX   = 6.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: SEVEN = 7.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: EIGHT = 8.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: NINE  = 9.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: TEN   = 10.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: HALF  = ONE/TWO
      REAL(KIND=DBLE), PARAMETER :: PI    = 3.14159265358979_DBLE
!
     CONTAINS 
!
    SUBROUTINE KEY_SORT(LIST,PERMUTATION)
!
!     This subroutine performs a key sort on the real LIST by the
!     heap sort method.  The returned permutation has the property that
!     LIST(PERMUTATION(I))<=LIST(PERMUTATION(I+1)) for all relevant I.
!     See: Nijenhuis A and Wilf HS (1978) "Combinatorial Algorithms for
!     Computers and Calculators, 2nd Ed.", Academic Press.
!
     IMPLICIT NONE
     INTEGER :: I,J,K,L,N,PSTAR
     INTEGER, DIMENSION(:) :: PERMUTATION
     REAL(KIND=DBLE), DIMENSION(:) :: LIST
!
!     Initialize the permutation.
!
     N = SIZE(LIST)
     DO I = 1,N
        PERMUTATION(I) = I
     END DO
     IF (N<=1) RETURN
!
!     Carry out the heap sort on the permutation key.
!
     L = 1+N/2
     K = N
     DO
        IF (L>1) THEN
           L = L-1
           PSTAR = PERMUTATION(L)
        ELSE
           PSTAR = PERMUTATION(K)
           PERMUTATION(K) = PERMUTATION(1)
           K = K-1
           IF (K==1) THEN
              PERMUTATION(1) = PSTAR
              RETURN
           END IF
        END IF
        I = L
        J = L+L
        DO WHILE (J<=K)
           IF (J<K) THEN
              IF (LIST(PERMUTATION(J))<LIST(PERMUTATION(J+1))) J = J+1
           END IF
           IF (LIST(PSTAR)<LIST(PERMUTATION(J))) THEN
              PERMUTATION(I) = PERMUTATION(J)
              I = J
              J = J+J
           ELSE
              J = K+1
           END IF
        END DO
        PERMUTATION(I) = PSTAR
     END DO
     END SUBROUTINE KEY_SORT
!
!
      SUBROUTINE FETCH_PREDICTOR_VALUES(X_INPUT,X,M,PEOPLE,PARAMETERS)
!
!     This subroutine loads the values of predictor M into the
!     vector X.
!
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,N,PEOPLE,PARAMETERS
      INTEGER(KIND=1) :: PHENOTYPE
      REAL(KIND=DBLE), DIMENSION(:) :: X
      REAL(KIND=DBLE), DIMENSION(PEOPLE,PARAMETERS) :: X_INPUT
!
!     Predictor number 1 is the grand mean.
!
      IF (M==1) THEN
         X = ONE
      ELSE
         X = REAL(X_INPUT(:,M),KIND=DBLE)
      END IF
      END SUBROUTINE FETCH_PREDICTOR_VALUES
!
!
!
!
    END MODULE COORDINATE_DESCENT
!
!


      SUBROUTINE LASSO_PENALIZED_L2_REGRESSION(X_INPUT,Y,LAMBDA,CASES,PREDICTORS,L2,R,OBJECTIVE,PENALTY,ESTIMATE)
!
!     New version written on 04/01/2007 by Ken Lange by cyclic coordinate descent.
!     This subroutine carries out penalized L2 regression with design
!     matrix X, dependent variable Y, and penalty constant LAMBDA.
!     Note that the rows of X correspond to cases and the columns
!     to predictors.  The first predictor is the constant 1.
!
      USE COORDINATE_DESCENT
      IMPLICIT NONE
      INTEGER :: CASES,I,ITERATION,J,K,L,M,N,PREDICTORS
      REAL(KIND=DBLE) :: CRITERION = TEN**(-5),EPSILON = TEN**(-8)
      REAL(KIND=DBLE) :: A,B,C,DL2,LAMBDA,L2,OBJECTIVE,PENALTY,NEW_OBJECTIVE
      REAL(KIND=DBLE) :: LEFT_L2,LEFT_OBJECTIVE,LEFT_PENALTY,LEFT_ROOT
      REAL(KIND=DBLE) :: RIGHT_L2,RIGHT_OBJECTIVE,RIGHT_PENALTY,RIGHT_ROOT
      REAL(KIND=DBLE), DIMENSION(CASES) :: Y
      REAL(KIND=DBLE), DIMENSION(PREDICTORS) :: ESTIMATE
      REAL(KIND=DBLE), DIMENSION(PREDICTORS,CASES) :: X_INPUT
      REAL(KIND=DBLE), DIMENSION(SIZE(X_INPUT,1)) :: SUM_X_SQUARES
      REAL(KIND=DBLE), DIMENSION(SIZE(Y)) :: LEFT_R,R,RIGHT_R
      REAL(KIND=DBLE), DIMENSION(SIZE(X_INPUT,2),SIZE(X_INPUT,1)) ::
      X = TRANSPOSE(X_INPUT)
!
!     Check that the number of cases is well defined.
!
!      IF (SIZE(Y)/=SIZE(X,1)) THEN
!        RETURN
!        PRINT*," THE NUMBER OF CASES IS NOT WELL DEFINED."
!      END IF
!
!     Check that the leftmost column of X consists of 1's.
!
!      IF (SUM(ABS(X(:,1)-ONE))>EPSILON) THEN
!         PRINT*," TOP ROW OF X SHOULD CONSIST OF 1'S"
!         RETURN
!      END IF
!
!     Initialize the number of cases M, the number of regression
!     coefficients N, and the sum of squares.
!
      M = SIZE(Y)
      N = SIZE(X,2)
      SUM_X_SQUARES = ZERO
!
!     Initialize the residual vector R and the penalty.
!
      R = Y
      IF (ABS(ESTIMATE(1))>ZERO) R = R-ESTIMATE(1)
      PENALTY = ZERO
      DO I = 2,N
         A = ESTIMATE(I)
         B = ABS(A)
         IF (B>ZERO) THEN
            R = R-A*X(:,I)
            PENALTY = PENALTY+B
         END IF
      END DO
!
!     Initialize the objective function and penalty.
!
      L2 = SUM(R**2)
      PENALTY = LAMBDA*PENALTY
      OBJECTIVE = HALF*L2+PENALTY
      !PRINT*," ITERATION = ",0," FUN = ",OBJECTIVE
!
!     Enter the main iteration loop.
!
      DO ITERATION = 1,1000
!
!     Update the intercept.
!
         A = ESTIMATE(1)
         ESTIMATE(1) = A+SUM(R)/REAL(M,KIND=DBLE)
         R = R+A-ESTIMATE(1)
!
!     Update the other regression coefficients.
!
         DO I = 2,N
            DL2 = -SUM(R*X(:,I))
            A = ESTIMATE(I)
            B = ABS(A)
            IF (B<EPSILON) THEN
               IF (DL2+LAMBDA>=ZERO.AND.-DL2+LAMBDA>=ZERO) CYCLE
            END IF
!
!     Find the root to the right of 0.
!
            IF (SUM_X_SQUARES(I)<=ZERO) SUM_X_SQUARES(I) = SUM(X(:,I)**2)
            RIGHT_ROOT = MAX(A-(DL2+LAMBDA)/SUM_X_SQUARES(I),ZERO)
            RIGHT_L2 = ZERO
            C = A-RIGHT_ROOT
            DO J = 1,M
               RIGHT_R(J) = R(J)+C*X(J,I)
               RIGHT_L2 = RIGHT_L2+RIGHT_R(J)**2
            END DO
            RIGHT_PENALTY = PENALTY+LAMBDA*(RIGHT_ROOT-B)
            RIGHT_OBJECTIVE = HALF*RIGHT_L2+RIGHT_PENALTY
!
!     Find the root to the left of 0.
!
            LEFT_ROOT = MIN(A-(DL2-LAMBDA)/SUM_X_SQUARES(I),ZERO)
            LEFT_L2 = ZERO
            C = A-LEFT_ROOT
            DO J = 1,M
               LEFT_R(J) = R(J)+C*X(J,I)
               LEFT_L2 = LEFT_L2+LEFT_R(J)**2
            END DO
            LEFT_PENALTY = PENALTY+LAMBDA*(ABS(LEFT_ROOT)-B)
            LEFT_OBJECTIVE = HALF*LEFT_L2+LEFT_PENALTY
!
!     Choose between the two roots.
!
            IF (RIGHT_OBJECTIVE<=LEFT_OBJECTIVE) THEN
               R = RIGHT_R
               ESTIMATE(I) = RIGHT_ROOT
               L2 = RIGHT_L2
               PENALTY = RIGHT_PENALTY
            ELSE
               R = LEFT_R
               ESTIMATE(I) = LEFT_ROOT
               L2 = LEFT_L2
               PENALTY = LEFT_PENALTY
            END IF
         END DO
!
!     Output the iteration number and value of the objective function.
!
         NEW_OBJECTIVE = HALF*L2+PENALTY
         !IF (ITERATION==1.OR.MOD(ITERATION,10)==0) THEN
            !PRINT*," ITERATION = ",ITERATION," FUN = ",NEW_OBJECTIVE
         !END IF
!
!     Check for a descent failure or convergence.  If neither occurs,
!     record the new value of the objective function.
!
         IF (NEW_OBJECTIVE>OBJECTIVE) THEN
            RETURN
            PRINT*," *** ERROR *** OBJECTIVE FUNCTION INCREASE"
            EXIT
         END IF
         IF (OBJECTIVE-NEW_OBJECTIVE<CRITERION) THEN
            RETURN
         ELSE
            OBJECTIVE = NEW_OBJECTIVE
         END IF
      END DO
      END SUBROUTINE LASSO_PENALIZED_L2_REGRESSION
!
!
!
    SUBROUTINE LASSO_PENALIZED_ESTIMATION(X_INPUT,Y,LAMBDA,PEOPLE,PARAMETERS, &
    OBJECTIVE,LOGLIKELIHOOD,R,ESTIMATE)
!
!     This subroutine carries out lasso penalized regression with dependent
!     variable Y and penalty constant LAMBDA.  Ordinary regression is
!     performed when REGRESSION_MODEL equals 1 and logistic regression when
!     REGRESSION_MODEL equals 2.  X represents one column of the design
!     matrix, whose rows correspond to people and whose columns correspond
!     to predictors.  The first predictor is the constant 1. Whenever the
!     number of active predictors exceeds a preset maximum, the subroutine
!     aborts.
!
    USE COORDINATE_DESCENT
    IMPLICIT NONE
    INTEGER :: I,ITERATION,J,K,L,M,MAX_PREDICTORS,PREDICTORS,PEOPLE,PARAMETERS
    INTEGER :: CARDINALITY,N
    LOGICAL :: DONE
    REAL(KIND=DBLE) :: CRITERION,LAMBDA,LOGLIKELIHOOD,L2,OBJECTIVE,PENALTY
    REAL(KIND=DBLE) :: BIG = FIVE*TEN,CONVERGENCE = TEN**(-5),EPSILON = TEN**(-8)
    REAL(KIND=DBLE) :: A,B,C,D,E,F,MU,NEW_OBJECTIVE,ROOT_CRITERION
!    INTEGER, OPTIONAL, DIMENSION(:) :: ACTIVE
    INTEGER, DIMENSION(:), ALLOCATABLE :: MEMBER
    INTEGER, DIMENSION(PARAMETERS) :: ACTIVETEMP
    REAL(KIND=DBLE), DIMENSION(PEOPLE,PARAMETERS) :: X_INPUT
    REAL(KIND=DBLE) :: INFORMATION,NEW_ROOT,PROB,ROOT,SCORE
    REAL(KIND=DBLE), DIMENSION(SIZE(X_INPUT,2)) :: ESTIMATE
    REAL(KIND=DBLE), DIMENSION(SIZE(X_INPUT,1)) :: X,INNER,R,Y
!
!     Initialize the current number of predictors and the estimates.
!
    PREDICTORS = 0
    ESTIMATE = ZERO
!
!     Initialize the residual vector R, the PENALTY, the OBJECTIVE function,
!     and an array INNER holding either the sum of squared predictors for
!     each parameter in ordinary regression or the inner product of each
!     predictor vector with the parameter vector in logistic regression,
!
    R = Y-HALF
    INNER = ZERO
    LOGLIKELIHOOD = -PEOPLE*LOG(TWO)
    PENALTY = ZERO
    OBJECTIVE = LOGLIKELIHOOD-PENALTY
!
!     Enter the main iteration loop.
!
    DO ITERATION = 1,1000
!
!     M is the index of the current active parameter.
!
    !M = 1
!
!     Update each regression coefficient in turn.
!
        DO I = 1,PARAMETERS
!
           ! CALL NEXT_SUBSET(MEMBER,CARDINALITY,N,DONE)
!
!     Skip inactive parameters.
!
   !         IF (PRESENT(ACTIVE)) THEN
   !             IF (I/=ACTIVE(M)) THEN
    !              !IF (DONE) CARDINALITY = CARDINALITY+1
   !                     CYCLE
   !             END IF
   !             M = M+1
    !       END IF
!
!     Fetch the values of the current predictor.
!
            CALL FETCH_PREDICTOR_VALUES(X_INPUT,X,I,PEOPLE,PARAMETERS)
!
!     Set the penalty for the current parameter.
!
            IF (I==1) THEN
                MU = ZERO
            ELSE
                MU = LAMBDA
            END IF
!
!     If a parameter is zero and both directional derivatives are
!     nonpositive, then cycle.
!
            A = ESTIMATE(I)
            B = ABS(A)
            IF (B<EPSILON) THEN
                SCORE = DOT_PRODUCT(R,X(1:PEOPLE))
            IF (ABS(SCORE)<=MU) CYCLE
            END IF
!
!     Optimize the objective function with respect to the current parameter.
!     The variable REGRESSION_MODEL is 1 for ordinary regression and 2 for
!     logistic regression.
!
            ROOT = A
            DO K = 1,10
!
!     Compute the score and information.
!
                SCORE = ZERO
                INFORMATION = ZERO
                C = ROOT-A
                DO J = 1,PEOPLE
                    D = X(J)
                    E = INNER(J)+C*D
                    IF (E<-BIG) THEN
                        PROB = ZERO
                    ELSE IF (E>BIG) THEN
                        PROB = ONE
                    ELSE
                        E = EXP(E)
                        PROB = E/(ONE+E)
                    END IF
                        SCORE = SCORE+(Y(J)-PROB)*D
                        INFORMATION = INFORMATION+PROB*(ONE-PROB)*D*D
                END DO
!
!     Apply Newton's method with directional safeguards.
!
                IF (ROOT>=EPSILON) THEN
                    NEW_ROOT = MAX(ROOT+(SCORE-MU)/INFORMATION,ZERO)
                ELSE IF (ROOT<=-EPSILON) THEN
                    NEW_ROOT = MIN(ROOT+(SCORE+MU)/INFORMATION,ZERO)
                ELSE
                    IF (SCORE>MU) THEN
                        NEW_ROOT = ROOT+(SCORE-MU)/INFORMATION
                    ELSE IF (SCORE<-MU) THEN
                        NEW_ROOT = ROOT+(SCORE+MU)/INFORMATION
                    ELSE
                        NEW_ROOT = ZERO
                    END IF
                END IF
!
!     Exit when the covergence criterion is satisfied.
!
                CRITERION = ABS(NEW_ROOT-ROOT)
                ROOT = NEW_ROOT
                IF (CRITERION<=EPSILON) EXIT
            END DO
!
!     Find the corresponding loglikelihood, penalty, and objective function.
!
            LOGLIKELIHOOD = ZERO
            C = ROOT-A
            DO J = 1,PEOPLE
                D = INNER(J)+C*X(J)
                INNER(J) = D
                IF (D<-BIG) THEN
                    PROB = ZERO
                    LOGLIKELIHOOD = LOGLIKELIHOOD+Y(J)*D
                ELSE IF (D>BIG) THEN
                    PROB = ONE
                    LOGLIKELIHOOD = LOGLIKELIHOOD+Y(J)*D-D
                ELSE
                    E = EXP(D)
                    F = ONE+E
                    PROB = E/F
                    LOGLIKELIHOOD = LOGLIKELIHOOD+Y(J)*D-LOG(F)
                END IF
                    R(J) = Y(J)-PROB
            END DO
            PENALTY = PENALTY+MU*(ABS(ROOT)-B)
            NEW_OBJECTIVE = LOGLIKELIHOOD-PENALTY
!
!     Record the estimate, and decrement the number of predictors if necessary.
!     If only active parameters are updated, then keep them from going inactive.
!
            ESTIMATE(I) = ROOT
            IF (ABS(ROOT)<EPSILON) THEN
                IF (B>=EPSILON) PREDICTORS = PREDICTORS-1
                ELSE
                    IF (B<EPSILON) PREDICTORS = PREDICTORS+1
            END IF
        END DO
!
!     Output the iteration number and value of the objective function.
!
        !IF (ITERATION==1.OR.MOD(ITERATION,10)==0) THEN
            !PRINT*," ITERATION = ",ITERATION," FUN = ",NEW_OBJECTIVE,LOGLIKELIHOOD,PENALTY
        !END IF
!
!     Check for a descent failure or an ascent failure.
!
        IF (NEW_OBJECTIVE<OBJECTIVE-EPSILON) THEN
            RETURN
            PRINT*," *** ERROR *** OBJECTIVE FUNCTION DECREASE"
            EXIT
        END IF
!
!     Check for convergence.
!
        CRITERION = ABS(NEW_OBJECTIVE-OBJECTIVE)
        OBJECTIVE = NEW_OBJECTIVE
        IF (CRITERION<CONVERGENCE) THEN
            !PRINT*," FINAL ITERATION = ",ITERATION," FUN = ",NEW_OBJECTIVE,LOGLIKELIHOOD,PENALTY
            EXIT
        END IF
    END DO
!
!     Deallocate the membership array and return.
!
    IF (ALLOCATED(MEMBER)) DEALLOCATE(MEMBER)
!
END SUBROUTINE LASSO_PENALIZED_ESTIMATION