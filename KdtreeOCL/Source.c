#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS


#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_SOURCE_SIZE     0x10000       
#define ALIGNMENT			4096    //serve per il zero copy
#define MAX_DIM 128
#define MAX_RAND 100
#define DIMTESTA 1000000
#define DIMTESTB 2000
#define DIST_MAX 99999999
#define GROUPSIZE 2

//PROVA RICCARDO

typedef struct node{
	float x[MAX_DIM];
	int idxB;
	int cd;
	int left;
	int right;
	int up;
	int idxT;

}node_idx;

typedef struct node_ris{
	float best_dist;
	int idx_tree;
}node_R;

struct kdtree_idx{
	int idx;
	struct node_idx *vett;
	struct node_idx *tree;
	int D;
	int S;
};

void scambiaS(node_idx *a,  node_idx *b){

	node_idx tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

void printNode_idx( node_idx node){

	printf("(");
	for (int i = 0; i < MAX_DIM; i++)
		printf("%f, ", node.x[i]);
	printf(")-> idxB:%d up:%d l:%d r:%d\n", node.idxB, node.up, node.left, node.right);
	
}



int quicksort(node_idx *v, int first, int last, int cd){

	int i, j;
	node_idx *pivot = ( node_idx*) malloc(sizeof(node_idx));



	if (first < last){
		i = first;
		j = last;
		*pivot = v[(first + last) / 2];


		do{
			while (v[i].x[cd] < pivot->x[cd]) i++;
			while (v[j].x[cd] > pivot->x[cd]) j--;
			if (i <= j){
				scambiaS(&v[i], &v[j]);
				i++;
				j--;
			}
		} while (i <= j);
		quicksort(v, first, j, cd);
		quicksort(v, i, last, cd);
	}

	return (first + last) / 2;
}




int findMedianMioOCL(node_idx* start, node_idx* end, int cd){

	int first, last, median;

	if (end <= start) return NULL;
	if (end == start + 1)
		return 0;

	//devo dare a quicksort start e end in int

	first = 0;
	last = (end - start) - 1;

	median = quicksort(start, first, last, cd);

	//printRoot(start);


	//struct kd_node_t *md = start + (end - start) / 2;

	//return md;
	return median;


}

node_idx* insert_node( node_idx* vett, node_idx* tree, int median, int cd, int idx){
	/*	cout << "inserito punto " << idx << " con idxvett " << median << " val ";
	for (int i = 0; i < S; i++){
	cout << vett[median].val[i];
	}cout << " figlio di "<<tree[idx].up<< endl;
	*/
	//tree[idx] = vett[median];
	tree[idx].idxB = vett[median].idxB;
	tree[idx].cd = cd;
	//tree[idx].cd = cd;
	for (int i = 0; i < MAX_DIM; i++){
		//tree[idx].x[i] = vett[median].x[i];
		tree[idx].x[i] = vett[median].x[i];
	}

	return tree;
}


node_idx* assignUp(node_idx* v, int idx, int i){

	v[i].up= idx;
	return v;
	
}

node_idx* assignSon(node_idx* v, int idx, int i, int flag){

	//if (idx == 0 && )
	
	//right
	if (flag == 1){
		v[i].right = idx;
		//v[i].left = -1;

	}
	else{
		v[i].left = idx;
		//v[i].right = -1;
	}

	return v;

}

/**
* Quicksort.
* @param a - The array to be sorted.
* @param first - The start of the sequence to be sorted.
* @param last - The end of the sequence to be sorted.
*/
void quickSortP(int first, int last, int depth, node_idx* vett)
{
	int axis = depth % MAX_DIM;
	int pivotElement;

	if (first < last)
	{
		pivotElement = pivot(first, last, depth, vett);
		quickSortP(first, pivotElement - 1, depth, vett);
		quickSortP(pivotElement + 1, last, depth, vett);
	}
}


int do_quicksort( int begin,  int end,  int depth, node_idx* vett){
	quickSortP(begin, end, depth, vett);
	return (begin + end) / 2;

}



/**
* Find and return the index of pivot element.
* @param a - The array.
* @param first - The start of the sequence.
* @param last - The end of the sequence.
* @return - the pivot element
*/
int pivot( int first, int last, int depth, node_idx* vett)
{
	int axis = depth % MAX_DIM;
	int  p = first;
	float pivotElement = vett[first].x[axis];

	for (int i = first + 1; i <= last; i++)
	{
		
		if (vett[i].x[axis] <= pivotElement)
		{
			p++;
			scambiaS(&vett[i], &vett[p]);
		}
	}

	scambiaS(&vett[p], &vett[first]);

	return p;
}




//def idx
int idx;

void  createTree( node_idx* buffer, node_idx** tree, int begin, int end, int idx_loc, int depth){

	int median;
	int cd = depth % MAX_DIM;
	//struct kd_node_t *n = (struct kd_node_t*) malloc(sizeof(struct kd_node_t));

	

	//selezione la cutting dimension corrente
	//int axis = cd % MAX_DIM;
	//struct kd_node_t *median = (struct kd_node_t*) malloc(sizeof(struct kd_node_t));

	
	//median = findMedianMioOCL(buffer, buffer + lenB, cd);
	median = do_quicksort(begin, end, depth, buffer);
	
	//stampa vettore ordinato
	/*
	for (int i = 0; i < DIMTESTB; i++){
		printf("(");
		for (int j = 0; j < MAX_DIM; j++)
			printf("%f, ", buffer[i].x[j]);
		printf(") ");
	}
	printf("\n");
		*/
		*tree = insert_node(buffer, *tree, median, cd, idx);
		depth++;
		
		idx_loc = idx;
		//left
		if ((median - 1 >= begin)) {
			int son = idx + 1;
			*tree = assignUp(*tree, idx_loc, son);
			//*tree[son].up = idx_loc;
			*tree = assignSon(*tree, son, idx_loc, 0);
			//*tree[idx_loc]->left = son;
			idx++;
			createTree(buffer, tree, begin, median - 1, idx, depth);
			//node->left = make_tree(begin, median - 1, depth);
			//node->left->up = node;
		}
		//right
		if ((median + 1 <= end)) {
			int son = idx + 1;
			//tree[son]->up = idx_loc;
			*tree = assignUp(*tree, idx_loc, son);
			//tree[idx_loc]->right = son;
			*tree = assignSon(*tree, son, idx_loc, 1);
			idx++;
			createTree(buffer, tree, median + 1, end, idx, depth);
			//node->right = make_tree(median + 1, end, depth = depth);
			//node->right->up = node;
		}
		

}



/*
struct kd_node_t{
	float x[MAX_DIM];
	int idx;
};






float
dist(struct kd_node_t *a, struct kd_node_t *b, int dim)
{
	float t, d = 0;
	while (dim--) {
		t = a->x[dim] - b->x[dim];
		d += t * t;
	}
	return d;
}


void swap(struct kd_node_t *x, struct kd_node_t *y) {
	float tmp[MAX_DIM];
	memcpy(tmp, x->x, sizeof(tmp));
	memcpy(x->x, y->x, sizeof(tmp));
	memcpy(y->x, tmp, sizeof(tmp));
}

//SCAMBIA NODE
void scambiaS(struct kd_node_t *a, struct kd_node_t *b){

	struct kd_node_t tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}


struct kd_node_t* assign(struct kd_node_t* buffer, struct kd_node_t* node, int idx){


	for (int i = 0; i < MAX_DIM; i++)
		buffer[idx].x[i] = node->x[i];

	buffer[idx].idx = idx;


	return buffer;


}




//QUICKSORT MYCREATE
void quicksort(struct kd_node_t *v, int first, int last, int cd){

	int i, j;
	struct kd_node_t *pivot = (struct kd_node_t*) malloc(sizeof(struct kd_node_t));



	if (first < last){
		i = first;
		j = last;
		*pivot = v[(first + last) / 2];


		do{
			while (v[i].x[cd] < pivot->x[cd]) i++;
			while (v[j].x[cd] > pivot->x[cd]) j--;
			if (i <= j){
				scambiaS(&v[i], &v[j]);
				i++;
				j--;
			}
		} while (i <= j);
		quicksort(v, first, j, cd);
		quicksort(v, i, last, cd);
	}
}


//PRINT ROOT (vettore)
void printRoot(struct kd_node_t* root, int dim){

	for (int i = 0; i < dim; i++){

		printf("(");
		for (int j = 0; j < MAX_DIM; j++){

			printf("%f, ", root[i].x[j]);

		}

		printf(") ");

	}
	printf("\n");
}









/* global variable, so sue me */
/*
int visited;
int count;
int position;
int idx;






//BUBBLESORT STRUCT NODE
void bubblesortT(struct kd_node_t *v, int dim, int cd){

	int i;

	while (dim > 1){

		for (i = 0; i < dim - 1; i++)
		if (v[i].x[cd] > v[i + 1].x[cd]){
			scambiaS(&v[i], &v[i + 1]);
			dim--;

		}

		dim--;
	}


}





*/


/*
//ASSIGNMENT NODE
struct kd_node_t* assign(struct kd_node_t *a, struct kd_node_t *b){

	for (int i = 0; i < MAX_DIM; i++)
		a->x[i] = b->x[i];

	a->idx = b->idx;
	
	return a;

}

*/


/*
//MYCREATE

//FIND_MEDIAN MYCREATE
struct kd_node_t* findMedianMio(struct kd_node_t* start, struct kd_node_t * end, int cd){

	int first, last;

	if (end <= start) return NULL;
	if (end == start + 1)
		return start;

	//devo dare a quicksort start e end in int

	first = 0;
	last = (end - start) - 1;

	quicksort(start, first, last, cd);

	//printRoot(start);


	struct kd_node_t *md = start + (end - start) / 2;

	return md;


}



struct kd_node_t* createTree(struct kd_node_t* bufferStart, struct kd_node_t** bufferFinished, int lenB, int cd){


	struct kd_node_t *n = (struct kd_node_t*) malloc(sizeof(struct kd_node_t));

	if (!lenB) return 0;

	//selezione la cutting dimension corrente
	//int axis = cd % MAX_DIM;
	//struct kd_node_t *median = (struct kd_node_t*) malloc(sizeof(struct kd_node_t));

	if ((n = findMedianMio(bufferStart, bufferStart + lenB, cd))) {
		n->idx = idx;
		//*bufferFinished[idx] = *n;
		*bufferFinished = assign(bufferFinished, n, idx);
		cd = (cd + 1) % MAX_DIM;

		//left
		idx = idx * 2 + 1;
		n = createTree(bufferStart, bufferFinished, n - bufferStart, cd);
		//se ritorna ritorniamo al valore precedente di idx
		idx = (idx - 1) / 2;
	
		
		
		//right
		idx = idx * 2 + 2;
		n = createTree(bufferStart, *bufferFinished, n - bufferStart, cd);
		//se ritorna ritorniamo al valore precedente di idx
		idx = (idx - 1) / 2;
		
	}
	return n;

}




*/



//PROVA SEARCH ITERATIVO

typedef struct point_query{
	float x[MAX_DIM];
	int best_idx;
	float best_dist;


}pq;

typedef struct stack_s{
	int idx;
	int left;
	short int count;
	struct stack_s *prec;



}stack;

stack* newStack( stack* s){

	stack *new_st = s;
	new_st->idx = -1;
	new_st->left = 0; //false
	new_st->prec = NULL;
	new_st->count = 0;

	return new_st;
}

pq* newPointQuery( pq* p){

	pq *new_pq = p;
	new_pq->best_idx = -1;
	new_pq->best_dist = 9999999;

	return new_pq;
}


stack *push_to_stack(stack *st, int idx, int b, stack *m, int count){

	//stack *new_st = (stack*)malloc(sizeof(stack));
	
	stack *new_st = &m[count];
	new_st = newStack( new_st);
	new_st->idx = idx;
	new_st->left = b;
	new_st->prec = st;
	//	cout << "inserito " << new_st->idx << endl;
	return new_st;
}

stack* pop_from_stack(stack *st){
	if ((st == NULL) || (st->prec == NULL)) return NULL;
	//stack *tmp = st->prec;
	//	if (st->prec == NULL)cout << "dddddddddddvdfvd\n";
	st = st->prec;
	//	if (st == NULL)cout << "vdfvd\n";
	//	cout << "tornato " << st->idx << endl;
	st->count++;
	return st;

}

void do_dist_it(const int p, pq* pointS, node_idx *tree, int DIMTREE){
	if ((p == -1) || (p >= DIMTREE)) return;
	float dist = 0;
	for (int i = 0; i < MAX_DIM; i++){
		//dist += pow((pq->val[i]-p->val[i]), 2);
		dist += (pointS->x[i] - tree[p].x[i]) * (pointS->x[i] - tree[p].x[i]);
	}
	//dist = sqrt(dist);		
	if (dist < pointS->best_dist){
		pointS->best_dist = dist;
		pointS->best_idx = tree[p].idxB;
	}
}

void search3_it(pq* pointS, int p, node_idx* tree, stack *m){ //p = 0
	int flag = 1;

	//stack *root = (stack*)malloc(sizeof(stack));
	//prova
	stack *root = &m[0];
	root = newStack(root);
	stack *st = root;
	int count = 0;

	while (1){
		//	printf("inside search 3_cpu\n");
		while ((p != -1) && (p < DIMTESTB)){ //da aggiustare questa perchè non è -1 ma 0 se non siamo in root
		//while ((p < DIMTESTB) && flag ){
			
			
			int axis = tree[p].cd;
			int left = 0;  //false
			do_dist_it(p, pointS, tree, DIMTESTB);
			if (pointS->x[axis] <= tree[p].x[axis]){
				left = 1;  //true
				count++;
				st = push_to_stack(st, p, left, m, count);
				int tmp = p;
				p = tree[tmp].left;
				//if (p == 0)
				//	flag = 0;
					
				//	cout << p << endl;
				//search3(pq, tmp);
			}
			else{
				count++;
				st = push_to_stack(st, p, left, m, count);
				int tmp = p;
				p = tree[tmp].right;
				
				//	cout << p << endl;
				//search3(pq, tmp);
			}
		}
		st = pop_from_stack(st);
		while ((st != NULL) && (st->count >= 2)){
			st = pop_from_stack(st);
		}
		if (st != NULL){
			p = st->idx;
			if (p != -1){
				int axis = tree[p].cd;
				if ((pointS->x[axis] - tree[p].x[axis]) * (pointS->x[axis] - tree[p].x[axis]) < pointS->best_dist){
					//if (abs(pq->val[axis] - p->val[axis])<pq->best_dist){
					//if (pow(pq->val[axis] - p->val[axis],2)<pq->best_dist){
					//if (sqrt(pow(pq->val[axis] - p->val[axis], 2))<pq->best_dist){
					//int tmp = left ? tree[p].right : tree[p].left;
					//search3(pq, tmp);
					p = st->left ? tree[p].right : tree[p].left;
				}
				else{ 
					p = -1;
				}
			}
		}
		else break;
	}
}

















//#define rand1() (rand() / (double)RAND_MAX)
//#define rand_pt(v) { v.x[0] = rand1(); v.x[1] = rand1(); v.x[2] = rand1(); }
#define rand1() (float)(rand() / (float)(RAND_MAX / MAX_RAND))
//#define rand_pt(v) { v.x[0] = rand1(); v.x[1] = rand1(); v.x[2] = rand1(); }

node_idx random(node_idx v, int dim){
	float tmp;

	for (int i = 0; i < dim; i++){
		tmp = rand1();

		v.x[i] = tmp;
	}

	return v;



}







int main(void)
{
	int i, j, k, size, dimA, dimTree, a;
	float best_dist;


	cl_ulong time_start, time_end;
	double total_time;

	double start = 0, end = 0, elapsed = 0;
	


	cl_event event;
	cl_mem a_mem_obj_pinned = NULL;
	cl_mem tree_mem_obj_pinned = NULL;
	cl_mem ris_mem_obj_pinned = NULL;
	cl_mem malloc_mem_obj = NULL;
	cl_mem pqMatr_mem_obj = NULL;

	cl_program program = NULL;

	dimA = DIMTESTA;
	dimTree = DIMTESTB;
	size = MAX_DIM;
	best_dist = 100000;
	int count = 0;

	const size_t global_size = DIMTESTA;

	





	//printf("Matrice A %dx%d   Matrice B %dx%d\n\n", DIMTESTA, MAX_DIM, DIMTESTB, MAX_DIM);

	//struct kd_node_t *foundG = (struct kd_node_t *) malloc(gDim)

	//struct kd_node *foundG = (struct kd_node *) malloc(2 * sizeof(struct kd_node*));

	//test BRUTEFORCE

	//node_idx *AB = (node_idx *)malloc(DIMTESTA * sizeof(node_idx));



	//USO PRIME
	//node_idx *BB = (node_idx *)malloc(DIMTESTB * sizeof(node_idx));


	



	//TEST MIO

	//MATRICE A come vettore di struct nodi
	printf("A: %d \n", DIMTESTA);
	//node_idx* A = (node_idx*) _aligned_malloc(DIMTESTA * sizeof(node_idx), ALIGNMENT);
	node_idx* A = (node_idx*) calloc(DIMTESTA, sizeof(node_idx));

	for (i = 0; i < DIMTESTA; i++) {
		//rand_pt(A[i]);
		A[i] = random(A[i], MAX_DIM);

		//AB[i] = A[i]; //inizializzo AB come A

		/*
		printf("{");
		for (j = 0; j < MAX_DIM; j++)
			printf("%f, ", A[i].x[j]);
		printf("}");
		*/
	}
	printf("\n");

	//A[0].x[0] = 3.190;
	//printRoot(A, DIMTESTA);
	//	printf("\n");



	//CREAZIONE KDTREE B
	//node_idx *B = (node_idx*) _aligned_malloc(DIMTESTB * sizeof(node_idx), ALIGNMENT);
	node_idx *B = (node_idx*)calloc(DIMTESTB, sizeof(node_idx));

	//ALBERO SU ARRAY
	//node_idx *TREE = (node_idx*) _aligned_malloc(DIMTESTB * sizeof(node_idx), ALIGNMENT);
	 node_idx *TREE = (node_idx*) calloc(DIMTESTB, sizeof(node_idx));

	//VETTORE RISULTATO
	//node_R *RIS = (node_R*) _aligned_malloc(DIMTESTA * sizeof(node_R), ALIGNMENT);
	node_R *RIS = (node_R*)calloc(DIMTESTA, sizeof(node_R));

	//prova senza malloc
	stack *m = (stack*) malloc(sizeof(stack) * DIMTESTB);

	//prova tracciaB

	node_idx *PRIME = (node_idx*)calloc(DIMTESTB, sizeof(node_idx));



	

	printf("KDTREE B: %d\n", DIMTESTB);
	for (i = 0; i < DIMTESTB; i++) {
		//	rand_pt(B[i]);
		B[i] = random(B[i], MAX_DIM);
		PRIME[i] = B[i];
		B[i].idxB = i;
		PRIME[i].idxB = i;
		//B[i].cd = 0;
		//B[i].left = -1;
		//B[i].right = -1;
		TREE[i].left = -1;
		TREE[i].right = -1;
		
		//B[i].up = 0;
		/*
		printf("idx:%d {", B[i].idxB);
		for (j = 0; j < MAX_DIM; j++)
			printf("%f, ", B[i].x[j]);
		printf("} ");
		
		*/
	}
	printf("\n\n");

	
	//CREAZIONI MATRICI C E RIS

	float **C = (float **)malloc(DIMTESTA * sizeof(float *));

	for (i = 0; i < DIMTESTA; ++i)
		C[i] = (float *)malloc(DIMTESTB * sizeof(float));

	int *RISB = (int *)malloc(DIMTESTA * sizeof(int));
	
	

	//CREAZIONE KDTREE_BUFFER

	/*

	idx = 0;

	//void  createTree(struct node_idx* buffer, struct node_idx** tree, int begin, int end, int idx_loc, int depth)
	createTree(B, &TREE, 0, DIMTESTB - 1, idx, 0);

	for (i = 0; i < DIMTESTB; i++){
		
		TREE[i].idxT = i;
		//printf("idxT:%d ", TREE[i].idxT);
		//printNode_idx(TREE[i]);
	}

	*/
	/*
	int u, v;
	printf("B:\n");
	for (u = 0; u < DIMTESTB; u++){

		printf("%f ", B[u].x[0]);
	}


	printf("\nTREE:\n");
	for (v = 0; v < DIMTESTB; v++){

		printf("%f ", TREE[v].x[0]);
	}
	printf("\n");


	printf("PRIME:\n");
	for (v = 0; v < DIMTESTB; v++){

		printf("%f ", PRIME[v].x[0]);
	}
	printf("\n");

	*/


	




	//PROVA BRUTE FORCE

	//bruteForceC(AB, PRIME, RISB, C);

	/*
	printf("\nRISB: \n");
	for (int y = 0; y < DIMTESTA; y++){

		printf("%d ", RISB[y]);
	}
	printf("\n");
	*/

	/*
	
	//PROVA RICERCA ITERATIVO -> FUNZIONA
	for (int a = 0; a < DIMTESTA; a++){

		pq* pointS = (pq*)malloc(sizeof(pq));
		pointS = newPointQuery(pointS);

		for (int l = 0; l < MAX_DIM; l++)
			pointS->x[l] = A[a].x[l];

		search3_it(pointS, 0, TREE, m);

		printf("searching node : (");

		for (int m = 0; m < MAX_DIM; m++)
			printf("%f, ", pointS->x[m]);
		printf(")\n");

		//best idx è la posizione del punto nella matrice B non l'idx dell'albero 
		printf("best node idxB and dist: %d - %f\n", pointS->best_idx, pointS->best_dist);

		RIS[a].idx_tree = pointS->best_idx;

		free(pointS);

	}

	

	*/
	


	//BISOGNA CREARE FUORI LA MATRICE DI POINT QUERY E PASSARLA AL KERNEL

	/*

	pq * pointQueryMatr = (pq*) malloc(sizeof(pq)* DIMTESTA);


	//inizializzo la matrice

	for (a = 0; a < DIMTESTA; a++){
		//pointQueryMatr[a] = newPointQuery(&pointQueryMatr[a]);
		pointQueryMatr[a].best_idx = -1;
		pointQueryMatr[a].best_dist = 9999999;
		for (int l = 0; l < MAX_DIM; l++)
			pointQueryMatr[a].x[l] = A[a].x[l];


	}

	*/
	/*
	printf("\n");

	
	for (int m = 0; m < DIMTESTA; m++){
		printf("searching node : { ");
		for (int n = 0; n < MAX_DIM; n++){
			printf("%f, ", pointQueryMatr[m].x[n]);
		}
		printf("}\n");
				
		}
			
		*/
	

	
	pq * pointQueryMatr = (pq*)malloc(sizeof(pq)* DIMTESTA);
	//INIZIO OPENCL

	start = clock();

	FILE *fp = fopen("Search_it.cl", "r");
	if (!fp) exit((printf("Failed to load kernel.\n"), 1));
	char* source_str = (char*)malloc(MAX_SOURCE_SIZE);
	size_t source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose(fp);


	// get platform and device information
	cl_platform_id platform_id = NULL;
	cl_device_id device_id = NULL;
	cl_context context = NULL;
	cl_command_queue command_queue = NULL;
	cl_int ret = clGetPlatformIDs(1, &platform_id, NULL);
	if (ret) exit((printf("clGetPlatformIDs failed\n"), 1));
	ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, NULL);
	if (ret) exit((printf("clGetDeviceIDs failed\n"), 1));


	// create an OpenCL context
	context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
	if (ret) exit((printf("clCreateContext failed\n"), 1));


	// create a command queu
	command_queue = clCreateCommandQueue(context, device_id, NULL, &ret);
	if (ret) exit((printf("clCreateCommandQueue failed\n"), 1));

	

	//PINNED BUFFER
	/*
	a_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, DIMTESTA * sizeof(node_idx), A, &ret);
	if (ret) exit((printf("clCreateBufferA failed\n"), 1));

	tree_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, DIMTESTB * sizeof(node_idx), TREE, &ret);
	if (ret) exit((printf("clCreateBufferB failed\n"), 1));

	ris_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, DIMTESTA * sizeof(node_R), 0, &ret);
	if (ret) exit((printf("clCreateBufferC failed\n"), 1));
	*/

	//NO PINNED BUFFER

	/*
	// create memory buffers on the device for each vector 
	a_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY, DIMTESTA * sizeof(node_idx), NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	
	//pointQueryMatr
	pqMatr_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(pq)* DIMTESTA, NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));


	tree_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY, DIMTESTB * sizeof(node_idx), NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	ris_mem_obj_pinned = clCreateBuffer(context, CL_MEM_WRITE_ONLY, DIMTESTA * sizeof(node_R), NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	malloc_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(struct stack_s) * DIMTESTA, NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	pq * pointQueryMatr = (pq*)malloc(sizeof(pq)* DIMTESTA);


	//inizializzo la matrice

	for (a = 0; a < DIMTESTA; a++){
		//pointQueryMatr[a] = newPointQuery(&pointQueryMatr[a]);
		pointQueryMatr[a].best_idx = -1;
		pointQueryMatr[a].best_dist = 9999999;
		for (int l = 0; l < MAX_DIM; l++)
			pointQueryMatr[a].x[l] = A[a].x[l];


	}

	// copy the lists pointQueryMatr and tree to their respective memory buffers
	ret = clEnqueueWriteBuffer(command_queue, pqMatr_mem_obj, CL_TRUE, 0,
		DIMTESTA * sizeof(pq), pointQueryMatr, 0, NULL, NULL);

	ret |= clEnqueueWriteBuffer(command_queue, tree_mem_obj_pinned, CL_TRUE, 0,
		DIMTESTB * sizeof(node_idx), TREE, 0, NULL, NULL);
	if (ret) exit((printf("clEnqueueWriteBuffer failed\n"), 1));
	*/

	// create a program from the kernel source
	program = clCreateProgramWithSource(context, 1,
		(const char **)&source_str,
		(const size_t *)&source_size,
		&ret);
	if (ret) exit((printf("clCreateProgramWithSource failed\n"), 1));

	// build the program
	ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
	if (ret) {
		size_t len;
		char buffer[4096];
		printf("clBuildProgram failed\n");
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		exit(1);
	}

	cl_kernel kernel = NULL;
	// create an OpenCL kernel
	kernel = clCreateKernel(program, "nearest", &ret);
	if (ret) exit((printf("clCreateKernel failed\n"), 1));



	//NO PINNED BUFFER

	//start = clock();

	//creazione albero
	createTree(B, &TREE, 0, DIMTESTB - 1, idx, 0);

	for (i = 0; i < DIMTESTB; i++){

		TREE[i].idxT = i;
		//printf("idxT:%d ", TREE[i].idxT);
		//printNode_idx(TREE[i]);
	}

	// create memory buffers on the device for each vector 
	a_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY, DIMTESTA * sizeof(node_idx), NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	
	//pointQueryMatr
	pqMatr_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(pq)* DIMTESTA, NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));


	tree_mem_obj_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY, DIMTESTB * sizeof(node_idx), NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	ris_mem_obj_pinned = clCreateBuffer(context, CL_MEM_WRITE_ONLY, DIMTESTA * sizeof(node_R), NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	malloc_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(struct stack_s) * DIMTESTA, NULL, &ret);
	if (ret) exit((printf("clCreateBuffer failed\n"), 1));

	


	//inizializzo la matrice

	for (a = 0; a < DIMTESTA; a++){
		//pointQueryMatr[a] = newPointQuery(&pointQueryMatr[a]);
		pointQueryMatr[a].best_idx = -1;
		pointQueryMatr[a].best_dist = 9999999;
		for (int l = 0; l < MAX_DIM; l++)
			pointQueryMatr[a].x[l] = A[a].x[l];


	}

	// copy the lists pointQueryMatr and tree to their respective memory buffers
	ret = clEnqueueWriteBuffer(command_queue, pqMatr_mem_obj, CL_TRUE, 0,
		DIMTESTA * sizeof(pq), pointQueryMatr, 0, NULL, NULL);

	ret |= clEnqueueWriteBuffer(command_queue, tree_mem_obj_pinned, CL_TRUE, 0,
		DIMTESTB * sizeof(node_idx), TREE, 0, NULL, NULL);
	if (ret) exit((printf("clEnqueueWriteBuffer failed\n"), 1));





	//settaggio parametri kernel

	ret |= clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&ris_mem_obj_pinned);
	ret |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&pqMatr_mem_obj);
	ret |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&tree_mem_obj_pinned);
	ret |= clSetKernelArg(kernel, 3, sizeof(int), (void *)&dimTree);
	ret |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&malloc_mem_obj);
	//ret |= clSetKernelArg(kernel, 5, sizeof(int), (void *)&dimA);

	//ret |= clSetKernelArg(kernel, 5, sizeof(float), (void *)&best_dist);
	//ret |= clSetKernelArg(kernel, 6, sizeof(int), (void *)&count);


	if (ret) exit((printf("clSetKernelArg failed\n"), 1));

	//esecuzione kernel
	ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_size, NULL,  //null come local size
    		0, NULL, &event);



	if (ret) exit((printf("clEnqueueNDRangeKernel failed\n"), 1));

	/*
	clWaitForEvents(1, &event);

	clFinish(command_queue);

	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

	total_time = time_end - time_start;
	printf("Tempo esecuzione Kernel  = %0.3f ms\n", (total_time / 1000000.0));
	
	*/
	//copia ris_pinned in RIS

	ret = clEnqueueReadBuffer(command_queue, ris_mem_obj_pinned, CL_TRUE, 0,
		DIMTESTA * sizeof(node_R), RIS, 0, NULL, NULL);
	if (ret) exit((printf("clEnqueueReadBuffer failed\n"), 1));

	end = clock();

	elapsed = (end - start) / (CLOCKS_PER_SEC / 1000);

	printf("tempo totale: %fms\n", elapsed);
	/*
	//verifico RIS
	
	
	
	int pos;
	int posT;
	printf("\n RIS: \n\n");
	for (int p = 0; p < DIMTESTA; p++){
		printf("ris[%d] idxB = %d -- node: {", p, RIS[p].idx_tree);

		pos = RIS[p].idx_tree;
		
		for (j = 0; j < MAX_DIM; j++)
			printf("%f, ", PRIME[pos].x[j]);
		printf("} dist: %f\n", RIS[p].best_dist);
		



	}
	
	
	

	printf("OK\n");

	
	
	

	// clean up


	//RELEASE KERNEL 
	ret = clReleaseKernel(kernel);



	//FREE A E B

	ret |= clReleaseMemObject(a_mem_obj_pinned);
	ret |= clReleaseMemObject(pqMatr_mem_obj);
	ret |= clReleaseMemObject(tree_mem_obj_pinned);
	ret |= clReleaseMemObject(malloc_mem_obj);
	ret |= clReleaseMemObject(ris_mem_obj_pinned);
	if (ret) exit((printf("resource1 release failed\n"), 1));

	free(A);
	free(B);
	free(TREE);
	free(RIS);
	free(pointQueryMatr);
	free(m);
	free(PRIME);


	//	ret |= clEnqueueUnmapMemObject(command_queue, c_mem_obj_pinned, cData, 0, NULL, NULL);
	//ret |= clEnqueueUnmapMemObject(command_queue, ris_mem_obj_pinned, RisData, 0, NULL, NULL);
	//	ret |= clReleaseMemObject(c_mem_obj_pinned); 

	//_aligned_free(C);
	//ret |= clReleaseMemObject(ris_mem_obj_pinned);
	ret |= clReleaseCommandQueue(command_queue);
	ret |= clReleaseContext(context);
	if (ret) exit((printf("resource release failed\n"), 1));

	/*
	endTOT = clock();

	elapsedTOT = (endTOT - startTOT) / (CLOCKS_PER_SEC / 1000);
	printf("\n");

	printf("Tempo esecuzione TOTALE = %f ms\n\n", elapsedTOT);

	*/
	//getchar();
	return 0;



	

}