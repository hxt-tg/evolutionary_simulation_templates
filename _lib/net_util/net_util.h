#ifndef NET_UTIL_H
#define NET_UTIL_H

#include <malloc.h>
#include <string.h>

// Debug
#include <stdio.h>

// Status
typedef char Status;
#define OK 0
#define FAILED 1
#define TRUE 1
#define FALSE 0

// Initial space
#define INCRE 5
#define INIT  10


// Player ID
typedef unsigned int PID;


// Shortcut: Default network Variable `net`.
#define NEI(n, a, bp) ((n)->p[(a)].nei[(bp)])
#define DEG(n, x) ((n)->p[(x)].count)

typedef struct SPlayer{
    PID *nei;    // Neighbors connection data
    int count;   // Count for neighbors
    size_t size; // Size of space
} *Player; 

typedef struct SNet {
    struct SPlayer *p;  //Struct Array, not pointer array.
    int n_players;	
    int n_links;
} *Net;

Net create_net(int n_player); //create a net with numbers of n_player
Status add_link(Net net, PID a, PID b); //
Status remove_link(Net net, PID a, PID b);
Status add_link_anyway(Net net, PID a, PID b); //
Status remove_link_anyway(Net net, PID a, PID b);
Status remove_link_by_position(Net net, PID a, int b_pos);
int degree(Net net, PID player);
int is_linked(Net net, PID a, PID b);
Status delete_net(Net net);

// For debug
void debug_print_net(Net net);
void debug_print_net_range(Net net, int start, int end);


#endif


