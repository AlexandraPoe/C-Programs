// ConsoleApplication4.cpp : This file contains the 'main' function. Program execution begins and ends there.
//                                   1st PDA LAboratory

/*#include <iostream>
#include <mpi.h>

#define MASTER 0

int main(int argc, char* argv[]) {
    int numprocs, procid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    if (numprocs != 2) {
        std::cerr << "This program requires exactly 2 processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int N = 0;
    if (procid == MASTER) {
        std::cout << "Enter the number of elements in the array (N): ";
        std::cin >> N;

        if (std::cin.fail() || N < 0 || N != static_cast<int>(N)) {
            std::cerr << "Number of elements must be a positive integer." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast N from master process to all other processes
    MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    // Allocate memory for the array
    int* array = new int[N];

    // Master process initializes and sends data
    if (procid == MASTER) {
        std::cout << "Enter " << N << " elements for the array:" << std::endl;
        for (int i = 0; i < N; ++i) {
            std::cin >> array[i] ;
            if (std::cin.fail() || array[i] < 0 || array[i] != static_cast<int>(N)) {
                std::cerr << "Elements must be positive integers." << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        int m;
        std::cout << "Enter the number to search (m): ";
        std::cin >> m;
        
        if (std::cin.fail() || m < 0 || m != static_cast<int>(m)) {
            std::cerr << "m must be a positive integer." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Send m and the array to the other process
        MPI_Send(&m, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(array, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    // Receiver process receives data and searches for m
    else {
        int m;
        MPI_Recv(&m, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(array, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Search for m in the array
        int position = -1;
        for (int i = 0; i < N; ++i) {
            if (array[i] == m) {
                position = i;
                break;
            }
        }

        if (position != -1) {
            std::cout << "Number " << m << " found at position " << position << " in the array." << std::endl;
        }
        else {
            std::cout << "Number " << m << " not found in the array." << std::endl;
        }
    }

    // Clean up
    delete[] array;
    MPI_Finalize();
    return 0;
}*/

                            //2nd PDA Laboratory: Collective Communication

/*1. Write a program that searches an element inside an array.
a. Use MPI_Broadcast for sending the array. If the element is found, print the maximum
position index. For computing the maximum position, you need to use MPI_Reduce.*/

/*
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>

#define MAXSIZE 100

int main(int argc, char** argv) {
    int myid, numprocs;
    std::vector<int> data(MAXSIZE);
    int myresult = 0, result = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (0 == myid) {
        std::ifstream fp("file.txt");
        if (!fp.is_open()) {
            std::cerr << "Can't open the input file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < MAXSIZE; i++) {
            fp >> data[i];
        }
        fp.close();
    }

    // Broadcast data
    MPI_Bcast(data.data(), MAXSIZE, MPI_INT, 0, MPI_COMM_WORLD);

    // Each process searches for a number
    int m;
    std::cout << "Enter the number to search: ";
    std::cin >> m;

    // Search for m in the array
    int position = -1;
   
    for (int i = 0; i < MAXSIZE; i++) {
        if (data[i] == m) {
            position = i;
        }
    }

    // Reduce to find the maximum position
    MPI_Reduce(&position, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // Output result
    if (position != -1) {
        if (myid == 0) {
            std::cout << "Number " << m << " found at position " << position << " in the array." << std::endl;
            std::cout << "The maximum position of " << m << " in the array is " << result << "." << std::endl;
        }
    }
    else {
        if (myid == 0) {
            std::cout << "Number " << m << " not found in the array." << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}*/

/*1. Write a program that searches an element inside an array.
b. Use scatter for sending the array. If the element is found many times, print all its positions.
Use MPI_Gather for sending back the positions.*/


/*#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>

#define MAXSIZE 100

int main(int argc, char** argv) {
    int myid, numprocs;
    std::vector<int> send_data(MAXSIZE);
    

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    std::vector<int> receive_data(MAXSIZE / numprocs);

    if (0 == myid) {
        std::ifstream fp("file1.txt");
        if (!fp.is_open()) {
            std::cerr << "Can't open the input file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < MAXSIZE; i++) {
            fp >> send_data[i];
        }
        fp.close();
    }

    // Scatter data
    MPI_Scatter(send_data.data(), MAXSIZE / numprocs, MPI_INT, receive_data.data(), MAXSIZE / numprocs, MPI_INT, 0, MPI_COMM_WORLD);

    // Each process searches for a number
    int m;
   std::cout << "Enter the number to search: ";
    std::cin >> m;

    // Search for m in the array
    int position = -1;
    std::vector<int> positions;
    for (int i = 0; i < MAXSIZE / numprocs; i++) {
        if (receive_data[i] == m) {
            positions.push_back(i + myid * MAXSIZE / numprocs);
        }
    }

    // Gather positions
    int positions_count = positions.size();
    std::vector<int> all_positions(numprocs * MAXSIZE / numprocs);
    MPI_Gather(&positions_count, 1, MPI_INT, all_positions.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(positions.data(), positions_count, MPI_INT, all_positions.data(), positions_count, MPI_INT, 0, MPI_COMM_WORLD);

    // Output result
    // Root process prints positions
    if (myid == 0) {
        std::cout << "Positions of number " << m << " in the array: ";
        for (int i = 0; i < all_positions.size(); i++) {
            if (all_positions[i] != 0) {
                std::cout << all_positions[i] << " ";
            }
        }
        std::cout << std::endl;
    }

    MPI_Finalize();
    return 0;
}*/

                        // 3rd PDA Laboratory: Derived Data Types

// The dataflow involves initializing MPI, defining a structured MPI datatype, 
// preparing and sending the data from rank 0 to all tasks, receiving the data in each task,
//  printing a sample of the received data in rank 0, and finalizing MPI
#include "mpi.h"
#include <iostream>
using namespace std;

#define NELEM 25

// Define the structure for particles
typedef struct {
    float x, y, z;
    float velocity;
    int n, type;
} Particle;

int main(int argc, char* argv[]) {
    int numtasks, rank, source = 0, dest, tag = 1, i;
    Particle p[NELEM], particles[NELEM];
    MPI_Datatype particletype, oldtypes[2];
    int blockcounts[2];
    MPI_Aint offsets[2], extent;
    MPI_Status stat;

    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    // Setup description of the 4 MPI_FLOAT fields x, y, z, velocity
    offsets[0] = 0;
    oldtypes[0] = MPI_FLOAT;
    blockcounts[0] = 4;

    // Setup description of the 2 MPI_INT fields n, type
    // Need to first figure offset by getting size of MPI_FLOAT
    MPI_Type_extent(MPI_FLOAT, &extent);
    offsets[1] = 4 * extent;
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 2;

    // Define structured type and commit it
    MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &particletype);
    MPI_Type_commit(&particletype);

    // Initialize the particle array and then send it to each task
    if (rank == 0) {
        for (i = 0; i < NELEM; i++) {
            particles[i].x = i * 1.0;
            particles[i].y = i * -1.0;
            particles[i].z = i * 1.0;
            particles[i].velocity = 0.25;
            particles[i].n = i;
            particles[i].type = i % 2;
        }
        for (i = 0; i < numtasks; i++)
            MPI_Send(particles, NELEM, particletype, i, tag, MPI_COMM_WORLD);
    }

    // Receive the particle array
    MPI_Recv(p, NELEM, particletype, source, tag, MPI_COMM_WORLD, &stat);

    // Print a sample of what was received
    if (rank == 0)
        for (i = 0; i < NELEM; i++)
            cout << "rank= " << rank << " " << p[i].x << " " << p[i].y << " " << p[i].z << " " << p[i].velocity << " " << p[i].n << " " << p[i].type << endl;

    // Free the MPI data type and finalize MPI
    MPI_Type_free(&particletype);
    MPI_Finalize();

    return 0;
}





// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
