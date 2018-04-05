iccg: iccg.c utils/vecmat.c
	gcc -o iccg.out iccg.c utils/vecmat.c

cg: cg.c utils/vecmat.c utils/mat_decomp.c
	gcc -o cg.out cg.c utils/vecmat.c utils/mat_decomp.c

st: stationary.c utils/vecmat.c
	gcc -o st.out stationary.c utils/vecmat.c

gauss: gaussian_elim.c
	gcc -o gauss.out gaussian_elim.c
