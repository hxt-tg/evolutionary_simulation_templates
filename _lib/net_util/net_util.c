#include "net_util.h"
#include <stdlib.h>

#define MALLOC(t) (t)Malloc(sizeof(struct S##t));

void *Malloc(size_t size){
    void *p = malloc(size);
    if (!p){
        printf("Malloc error.");
        exit(0);
    }
    return p;
}

void *Realloc(void *raw, size_t size){
    void *p = realloc(raw, size);
    if (!p){
        printf("Realloc error.");
        exit(0);
    }
    return p;
}

Net create_net(int n_players){
    Net net = MALLOC(Net);
    net->n_players= n_players;
    net->n_links = 0;
    net->p = (Player ) Malloc(n_players * sizeof(struct SPlayer));
    int i;
    for (i = 0; i < n_players; ++i) {
        net->p[i].size = INIT * sizeof(PID);
        net->p[i].count = 0;
        net->p[i].nei = (PID *) Malloc(net->p[i].size);
    }
    return net;
}

int is_linked(Net net, PID a, PID b){
    PID *i = net->p[a].nei;
    for (; i < net->p[a].nei + net->p[a].count; i++)
        if (b == *i)
        	return TRUE;
    return FALSE;
}

void _add_nei(Player p, PID id){
    if (p->size <= p->count * sizeof(PID)){
        p->nei = (PID*)Realloc(p->nei, p->size + INCRE * sizeof(PID));
        p->size += INCRE*sizeof(PID);
    }
    p->nei[p->count++] = id;
}

void _remove_nei(Player p, PID id){
    if (p->count == 0) return ;
    PID *i = p->nei;
    for (; i < p->nei + p->count; i++){
        if (id == *i){
            *i = *(p->nei + --p->count);
            break;
        }
    }
}

void _swap(PID *a, PID *b){
    *a = *a ^ *b;
    *b = *a ^ *b;
    *a = *a ^ *b;
}

Status add_link(Net net, PID a, PID b){
	if ((int)a >= net->n_players || (int)b >= net->n_players || a == b)
        return FAILED;
	if (!is_linked(net, a, b))
		return add_link_anyway(net, a, b);
    return FAILED;
}

Status add_link_anyway(Net net, PID a, PID b){
    // if (a > b) _swap(&a, &b);
    _add_nei(net->p + a, b);
    _add_nei(net->p + b, a);
    net->n_links ++;
    return OK;
}

Status remove_link(Net net, PID a, PID b){
	if ((int)a >= net->n_players || (int)b >= net->n_players || a == b)
        return FAILED;
    if (is_linked(net, a, b))
    	return remove_link_anyway(net, a, b);
    return FAILED;
}

Status remove_link_anyway(Net net, PID a, PID b){
    _remove_nei(net->p+a, b);
    _remove_nei(net->p+b, a);
    net->n_links --;
    return OK;
}

Status remove_link_by_position(Net net, PID a, int bp){
    Player pa = net->p + a;
    if ((int)a >= net->n_players || bp < 0 || bp >= net->p[a].count)
        return FAILED;
    pa->nei[bp] = pa->nei[pa->count-1];
    pa->count--;
    net->n_links --;
    return OK;
}

int degree(Net net, PID player){
    return net->p[player].count;
}

Status delete_net(Net net){
    if (!net->p) return FAILED;
    int i;
    for (i = 0; i < net->n_players; i++)
        free(net->p[i].nei);
    free(net->p);
    net->p = NULL;
    return OK;
}

// For debug
void debug_print_net(Net net){
	int i, j;
	for (i = 0; i < net->n_players; i++){
		printf("%d: ", i);
		for (j = 0; j < DEG(net, i)-1; j++)
			printf("%d, ", NEI(net, i, j));
		if (DEG(net, i) > 0)
			printf("%d", NEI(net, i, net->p[i].count-1));
		putchar(10);
	}
}

void debug_print_net_range(Net net, int start, int end){
	int i, j;
	for (i = start; i < (end > net->n_players ? net->n_players : end); i++){
		printf("%d: ", i);
		for (j = 0; j < DEG(net, i)-1; j++)
			printf("%d, ", NEI(net, i, j));
		if (DEG(net, i) > 0)
			printf("%d", NEI(net, i, net->p[i].count-1));
		putchar(10);
	}
}

