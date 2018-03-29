iccg: iccg.c utils/vecmat.c
	gcc -o iccg.out iccg.c utils/vecmat.c

cg: cg.c utils/vecmat.c
	gcc -o cg.out cg.c utils/vecmat.c
