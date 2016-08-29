#define MAX_DIM 1                        


typedef struct node{
	float x[MAX_DIM];
	int idx;
	int cd;
	int left;
	int right;
	int up;

}node_idx;

typedef struct node_ris{
	float x[MAX_DIM];
	int idx_tree;
}node_R;


int isEqual(node_idx a, node_idx b){ 
	

	for(int i = 0; i < MAX_DIM; i++){
		if(a.x[i] != b.x[i])
			return 0;
	}

	return 1;
	 
}



__kernel void nearest(__global node_R* RIS, __global node_idx* A, __global node_idx* TREE, int DIMTREE, int DIM_MAX, int DIST_BEST, int count){


	
 int idx = get_global_id(0); 
 
 int depth = 0;
 float d, dx, dx2, t;
 int k = 0;
 int cd = 0;
 float best = DIST_BEST;
 int i,j;
 int flag =  1;
 int countRoot = count;

 printf("%d\n", idx);
 
 while(flag){


		
	
	if( isEqual(TREE[k], TREE[0]) ){ 

		countRoot++;
		printf("%d\n", countRoot);
		
	
	 }

	 if(countRoot > 2){ 
		flag = 0;
		printf("fine\n");
	}
		else{

	 if(isEqual(TREE[k], TREE[0]) && countRoot >= 2 && TREE[0].right != -1){
		printf("->destra\n");
		k = TREE[0].right;
		}

	for(i = 0; i< DIM_MAX; i++){ 
		
		printf("%f - %f\n",A[idx].x[i], TREE[k].x[i]);
		t = A[idx].x[i] - TREE[k].x[i];
		d += t * t;
	}

	dx = TREE[k].x[cd] - A[idx].x[cd];
	dx2 = dx * dx;

	if(d < best ){ 
		printf("%f < best %f\n", d, best);

		for(j = 0; j< DIM_MAX; j++)
			RIS[idx].x[j] = TREE[k].x[i];

		RIS[idx].idx_tree = k;
		best = d;
		d = 0;
		printf("best : (%f)\n", TREE[k].x[0]); 
		
	}

	printf("dx:%f\n", dx);
	printf("dx2:%f\n", dx2);
	printf("best_d:%f\n", best);
	//albero sinistro
	if( dx > 0 && TREE[k].left != -1){ 
		printf("->alb sx\n");
	
		k = TREE[k].left;
		cd = TREE[k].cd;
		d = 0;
	//	depth++;
	 }
	 else if( dx <= 0 && TREE[k].right != -1){
		printf("->alb dx\n");
		k = TREE[k].right;
		cd = TREE[k].cd;
		d = 0;
	//	depth++;
	  }
	  else{
		 printf("up\n");
		 k = TREE[k].up; 
		 d = 0;
		 
		 }

	  //if(dx2 >= best)
		//k = TREE[k].up; 

		 
		 }
		 /*
		 if( countRoot + 1 > 2 ){
			flag = 0;
			printf("fine\n");
		}
		*/
		
			
		

		
	}
	

	
	
}


 
  



 