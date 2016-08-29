
//IMPLEMENTARE ZERO COPY

#define MAX_DIM 128
#define DIMTESTA 1000000 // == DIMTESTA
/*
typedef struct node{
	float x[MAX_DIM];
	int idx;
	int cd;
	int left;
	int right;
	int up;

}node_idx;

*/

struct node{
	float x[MAX_DIM];
	int idxB;
	int cd;
	int left;
	int right;
	int up;
	int idxT;

};

/*
typedef struct node_ris{
	float x[MAX_DIM];
	int idx_tree;
}node_R;

*/

struct node_ris{
	float best_dist;
	int idx_tree;
};

//PROVA SEARCH ITERATIVO
/*
typedef struct point_query{
	float x[MAX_DIM];
	int best_idx;
	float best_dist;


}pq;

*/

struct point_query{
	float x[MAX_DIM];
	int best_idx;
	float best_dist;


};

 struct stack_s{
	int idx;
	int left;
	short int count;
	__global struct stack_s *prec;



};

__global struct stack_s* newStack( __global struct stack_s* s){

	__global struct stack_s *new_st = s;
	new_st->idx = -1;
	new_st->left = 0; //false
	new_st->prec = NULL;
	new_st->count = 0;

	return new_st;
}
/*
struct point_query* newPointQuery( struct point_query* p){

	struct point_query *new_pq = p;
	new_pq->best_idx = -1;
	new_pq->best_dist = 9999999;

	return new_pq;
}

*/


__global struct stack_s *push_to_stack(__global struct stack_s *st, int idx, int b, __global struct stack_s* m, int count){

	//struct stack_s *new_st = (struct stack_s*)malloc(sizeof(struct stack_s));
	__global struct stack_s *new_st = &m[count];
	new_st = newStack( new_st);
	new_st->idx = idx;
	new_st->left = b;
	new_st->prec = st;
	//	cout << "inserito " << new_st->idx << endl;
	return new_st;
}

__global struct stack_s* pop_from_stack(__global struct stack_s *st){
	if ((st == NULL) || (st->prec == NULL)) return NULL;
	//stack *tmp = st->prec;
	//	if (st->prec == NULL)cout << "dddddddddddvdfvd\n";
	st = st->prec;
	//	if (st == NULL)cout << "vdfvd\n";
	//	cout << "tornato " << st->idx << endl;
	st->count++;
	return st;

}

//pointS era un puntatore
void do_dist_it(const int p, struct point_query *pointS,__global struct node *tree, int DIMTREE){
	if ((p == -1) || (p >= DIMTREE)) return;
	float dist = 0;
	for (int i = 0; i < MAX_DIM; i++){
		//dist += pow((pq->val[i]-p->val[i]), 2);
		dist += (pointS->x[i] - tree[p].x[i]) * (pointS->x[i] - tree[p].x[i]);  //->
	}
	//dist = sqrt(dist);		
	if (dist < pointS->best_dist){
		pointS->best_dist = dist;
		pointS->best_idx = tree[p].idxB;
	}
}

//era un puntatore pointS

void search3_it(struct point_query *pointS, int p,__global struct node* tree, int DIMTESTB, __global struct stack_s *m){ //p = 0
	int flag = 1;

	//struct stack_s *root = (struct stack_s*)malloc(sizeof(struct stack_s));
	__global struct stack_s *root = &m[0];
	root = newStack(root);
	__global struct stack_s *st = root;
	int count = 0;

	while (1){
		//	printf("inside search 3_cpu\n");
		while ((p != -1) && (p < DIMTESTB)){ //da aggiustare questa perchè non è -1 ma 0 se non siamo in root
		//while ((p < DIMTESTB) && flag ){
			
			
			int axis = tree[p].cd;
			int left = 0;  //false
			do_dist_it(p, pointS, tree, DIMTESTB);   //->da mutare in kernel
			if (pointS->x[axis] <= tree[p].x[axis]){
			//if (pointS.x[axis] <= tree[p].x[axis]){
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
				//if ((pointS.x[axis] - tree[p].x[axis]) * (pointS.x[axis] - tree[p].x[axis]) < pointS.best_dist){
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

__kernel void nearest(__global struct node_ris* RIS, __global struct point_query *PQ, __global struct node* TREE, int DIMTREE, __global struct stack_s* m){

 int idx = get_global_id(0); 
 int p = 0;

 struct point_query pq = PQ[idx];

 if(idx < DIMTESTA){ 
//printf("searching node pq %f %f %d \n", PQ[idx].x[0], PQ[idx].best_dist, PQ[idx].best_idx);
// struct point_query pq[DIMPQMATR];
//struct point_query *pq = &PQ[idx];

//printf("node2 pq %f %f %d \n", pq.x[0], pq.best_dist, pq.best_idx);
 //if(idx >= DIMPQMATR)
//	return;
// pq[idx] = PQ[idx]; 
 //printf("node pq %f %f %d \n", pq.x[0], pq.best_dist, pq.best_idx);
	search3_it(&pq, p, TREE, DIMTREE, m);
	}

	
//printf("searching node %f found %f %d \n", pq.x[0], pq.best_dist, pq.best_idx);
 RIS[idx].best_dist = pq.best_dist;
 RIS[idx].idx_tree = pq.best_idx;

//printf("searching node (%f)\nfound node idxB: %d - dist: %f\n",pq.x[0], RIS[idx].idx_tree, RIS[idx].best_dist);
//printf("node %f\n", TREE[RIS[idx].idx_tree].x[0]);

// printf("OK\n");
 


	

 






}


