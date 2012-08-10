
		case 25: //  Pressure vessel (penalty method)		
			pb.SS.D=4;
/*
		pb.SS.min[0] = 1.125; pb.SS.max[0] = 12.5;
		pb.SS.q.q[0] = 0.0625;
		pb.SS.min[1] = 0.625; pb.SS.max[1] = 12.5; 
		pb.SS.q.q[1] = 0.0625;
		pb.SS.min[2] = 0.00000001; pb.SS.max[2] = 240; 
		pb.SS.q.q[2] = 0;
		pb.SS.min[3] = 0.00000001; pb.SS.max[3] = 240; 
		pb.SS.q.q[3] = 0;
*/

		pb.SS.min[0] = 0.0625; pb.SS.max[0] = 99;
		pb.SS.q.q[0] = 0.0625;
		pb.SS.min[1] = 0.0625; pb.SS.max[1] = 99; 
		pb.SS.q.q[1] = 0.0625; 
		pb.SS.quantisation=1;
		
		pb.SS.min[2] = 10; pb.SS.max[2] = 200; 
		pb.SS.q.q[2] = 0;
		pb.SS.min[3] = 10; pb.SS.max[3] = 200; 
		pb.SS.q.q[3] = 0;


		pb.evalMax = 200000 ; 
		pb.epsilon = 0; //0.00001; //0.0000000001;			pb.objective =0; // 6059.714335 
		break;