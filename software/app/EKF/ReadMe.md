 EKF - Extended Kalman Filter for the Estimation Of Attitude of a Rigid Body

    Author: Paulo Henrique Vancin   

    The intended purpose of the implementation of this algorithm is to test
    capability of the riscv to run this kind of application. For now this
    code only runs offline getting the sensors data from a vector previously
    constructed. Later will be used for the control of a quadcopter.

    Inputs: - Accelerometer Data (3 components: X, Y and Z);
            - Gyroscope Data (3 components: X, Y and Z);
            - Magnetometer Data (3 components: X, Y and Z);
            - Time Stamp (Not Implemented Yet);
            - Barometer (Not Implemented Yet - Will not run in the EKF - only for Sync);


    Outputs: - Quaternions (4 components - later will be the euler angles);
             - Time Stamp (Not Implemented Yet);  
             - Barometer (Not Implemented Yet - Will not run in the EKF - only for Sync);


    Includes Needed: - "hf-risc.h" (Modded from https://github.com/sjohann81/hf-risc/tree/master/software);
                     - "fixed.h" (Modded from https://github.com/sjohann81/hf-risc/tree/master/software);
                     - "matrix.h";
                     - "quaternions.h";
                     - "setType.h";
                     - "EKF.h";

    HOW TO USE:

    - Set the initial values:
          - Ts: Period between sensors read (in seconds);
          - Matrix Q = Covariance Matrix of the Created Model (7X7 Matrix);
          - Matrix R = Covariance Matrix of the Sensors Read (4X4 Matrix);
          - Matrix P = Error Covariance Matrix (Usually a 7X7 Identity Matrix);
          - Set the initial value of the quaternions calling imu2q() with the
          first read from sensor;
          - Call the EKF() and get the first matrix for P and the States;

    - Get the sensors data for each sample;
    - Run the EKF for that sample;


    The "setType.h" include, grants the ability to choose the number format for
    the simulation (Fixed Point, Double or FLOAT). Put the parameter "-DT"Type""
    in the makefile to set the type.

    Ex: - Fixed Point: -DTFIXED;
        - Double: -DTDOUBLE;
        - Float: - DTFLOAT;


    Makefile Config:

    EKF: crt
    		$(GCC_$(ARCH)) -o fixed.o lib/fixed.c
    		$(GCC_$(ARCH)) -o matrix.o -c "folder"/matrix.c -w -s -g -DT"Type"
    		$(GCC_$(ARCH)) -o quaternions.o -c "folder"/quaternions.c -w -s -g -DT"Type"
    		$(GCC_$(ARCH)) -o EKF.o -c "folder"/EKF.c -w -s -g -DT"Type"
    		$(GCC_$(ARCH)) -o main.o -c "folder"/main.c -w -s -g -DT"Type"
    		@$(MAKE) --no-print-directory link

    BUGS: - Double not working due to an error found on "math.h", when a negative
          number is multiplied by zero, the resulted value is different than zero;
          - With Fixed Point, is needed to use covariance matrices with similar
          order of magnitude, otherwise the computation of the inverse matrix
          will result in a larger number than the fixed library is capable of
          represent;  
