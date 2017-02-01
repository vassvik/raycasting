#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Dist {
	char dx, dy, dz;
	unsigned char dist;

		bool operator<(const Dist &rhs) {
		return this->dist < rhs.dist;
	}
	Dist& operator=(const Dist &rhs) {
		this->dx = rhs.dx;
		this->dy = rhs.dy;
		this->dz = rhs.dz;
		this->dist = rhs.dist;

		return *this;
	}
};

#define MIN(a,b)   ((a) < (b) ? (a) : (b))
#define MAX(a,b)   ((a) > (b) ? (a) : (b))

const int Nx = 512, Ny = Nx, Nz = Nx;
unsigned char*data;

Dist *distances;

void imsort(Dist* xs, int l, int u);

void swap(Dist* xs, int i, int j) {
    Dist tmp = xs[i]; xs[i] = xs[j]; xs[j] = tmp;
}

// http://stackoverflow.com/questions/2571049/how-to-sort-in-place-using-the-merge-sort-algorithm
void wmerge(Dist* xs, int i, int m, int j, int n, int w) {
    while (i < m && j < n)
        swap(xs, w++, xs[i] < xs[j] ? i++ : j++);
    while (i < m)
        swap(xs, w++, i++);
    while (j < n)
        swap(xs, w++, j++);
} 

/* 
 * sort xs[l, u), and put result to working area w. 
 * constraint, len(w) == u - l
 */
void wsort(Dist* xs, int l, int u, int w) {
    int m;
    if (u - l > 1) {
        m = l + (u - l) / 2;
        imsort(xs, l, m);
        imsort(xs, m, u);
        wmerge(xs, l, m, m, u, w);
    }
    else
        while (l < u)
            swap(xs, l++, w++);
}

void imsort(Dist* xs, int l, int u) {
    int m, n, w;
    if (u - l > 1) {
        m = l + (u - l) / 2;
        w = l + u - m;
        wsort(xs, l, m, w); /* the last half contains sorted elements */

        while (w - l > 2) {
            n = w;
            w = l + (n - l + 1) / 2;
            wsort(xs, w, n, l);  /* the first half of the previous working area contains sorted elements */
            wmerge(xs, l, l + n - w, n, u, w);
        }
        for (n = w; n > l; --n) /*switch to insertion sort*/
            for (m = n; m < u && xs[m] < xs[m-1]; ++m)
                swap(xs, m, m - 1);
    }
}


void distanceEstimate(int num_search) {


	double avgK = 0.0;
	double avgDist = 0.0;
	double maxDist = -1.0;

	int count = 0;

	for (int k = 0; k < Nx*Ny*Nz; k++) {
		int xy = k % (Nx*Ny);
		int z = k / (Nx*Ny);
		int x = xy % Nx;
		int y = xy / Nx;


		if (data[k] == 0) 
			continue;

		double min_dist = 100000000;
		for (int K = 0; K < num_search; K++) {
			int X = x+distances[K].dx;
			int Y = y+distances[K].dy;
			int Z = z+distances[K].dz;
			// periodic boundaries
			if (X < 0)   X += Nx;
			if (X >= Nx) X -= Nx;
			if (Y < 0)   Y += Ny;
			if (Y >= Ny) Y -= Ny;
			if (Z < 0)   Z += Nz;
			if (Z >= Nz) Z -= Nz;
			int id = Z*Nx*Ny+Y*Nx+X;

			if (id == k)
				continue;

			if (data[id] == 0) {
				min_dist = distances[K].dist;
				avgK += K;
				count++;
				break;
			}
		}
		// KNOW maxDIST == 30 for segmented_castle_512.ubc and periodic boundaries
		data[k] = MIN(31, int(min_dist));
		avgDist += min_dist;
		maxDist = MAX(maxDist, min_dist);

		if ((k+1) % (20000) == 0) {
			printf("\r%d/%d (%.2f%%). Avg num searched = %f, avg dist = %.2f, max dist = %.2f      ",k+1, Nx*Ny*Nz, 100*(k+1.0)/(Nx*Ny*Nz), avgK/count, avgDist/count, maxDist); fflush(stdout);
		}
			
	}
	printf("\n"); fflush(stdout);
}

int main() {
	char filename[] = "../data/segmented_castle_512.ubc";
	int size = Nx*Ny*Nz;

	// Assuming file size is known
	printf("Reading file %s\n", filename); fflush(stdout);
	FILE *fp = fopen(filename, "rb");
    data = new unsigned char[size]; 
    fread(data, sizeof(unsigned char), size, fp);
    fclose(fp);

    printf("Pre-calculating distances\n"); fflush(stdout);
    // initialize distance calculations
    distances = new Dist[Nx*Ny*Nz];
    int k = 0;
    for (int z = 0; z < Nz/4; z++) {
		for (int y = 0; y < Ny/4; y++) {
				for (int x = 0; x < Nx/4; x++) {
					int dz = z;
					int dx = x;
					int dy = y;

					if (dx > Nx/8.0) dx -= Nx/4;
					if (dy > Ny/8.0) dy -= Ny/4;
					if (dz > Nz/8.0) dz -= Nz/4;

					// distance relative to (0,0,0)
					distances[k].dx = dx;
				distances[k].dy = dy;
				distances[k].dz = dz;
				distances[k].dist = (unsigned short)sqrt(dx*dx+dy*dy+dz*dz);
				k++;
			}
		}
	}

	// sort distances
	printf("Sorting relative distances (~3 sec)\n"); fflush(stdout);
	imsort(distances, 0, k);

	printf("Calculating voxel distances (~16 min)\n"); fflush(stdout);
	distanceEstimate(k);

	fp = fopen("../data/distance_segmented_castle_512.ubc", "wb");
    fwrite(data, sizeof(unsigned char), size, fp);
    fclose(fp);

    free(data);


	return 0;
}