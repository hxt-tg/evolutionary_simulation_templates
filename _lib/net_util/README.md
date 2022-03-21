# Net Utils

A simple implementation of complex network with dynamic neighbors by variable length of arrays.
Add/Remove links supported.
<span style="color: red">(WARNING: This version is unstable and with some unknown bugs.)</span>

## Author
- @hxt-tg:
	Email: hxt.taoge@gmail.com


## Shortcuts

- NEI(n, a, bp): The No.bp neighbor of node a in network n.
- DEG(n, x): The Degree of node x in network n.

## Status code

For some reason, this status code is UNIX-like, which means function will return 0 if all things are OK.
If some errors happen, a value greater than 0 will be returned.


## Example (Project-based required)

```cpp
#include <stdio.h>
#include "net_util.h"
#define SIZE 1000

int main(){
	Net net = create_net(SIZE);		// Create a net with 1000 nodes
	int i;
    for (i = 0; i < SIZE; i++)
        if (add_link(net, i, (i+1)%SIZE) > 0 || add_link(net, i, (i+2)%SIZE) > 0)
        	printf("Error occurred when producing network.\n");

    printf("Degree of node 100: %d\n", degree(net, 100));
    for (i = 0; i < DEG(net, 100); i++)
        printf("  The No.%d neighbor of node 100: %d\n", i, NEI(net, 100, i));
    printf("Are node 199 and 200 linked? %s.\n", is_linked(net, 199, 200) ? "Yes":"No");
    remove_link(net, 199, 200);
    printf("Are node 199 and 200 linked? %s.\n", is_linked(net, 199, 200) ? "Yes":"No");
    delete_net(net);
}
```