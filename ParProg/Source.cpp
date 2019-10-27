#include <mpi.h>
#include <iostream>
#include <queue>

using namespace std;

#define MANAGER 0
#define PRODUCER 2 
#define CONSUMER 1

struct info
{
	int rank;
	int whatDoYouNeed;
	int res;
};

const int PUT_RESOURCE = 1;
const int GET_RESOURCE = 2;
const int EXITP = 3;
const int EXITC = 4;
int STOP = -1;



class Consumer {

private:
	int resources_to_consume; 
	std::queue<int> resources; 
	info i1;
private:
	void RequestResource() { 

		MPI_Send(&i1, 3, MPI_INT, MANAGER, 0, MPI_COMM_WORLD); 
	}
	void RecieveResource() { 
		int resource = 0;
		MPI_Status status;
		MPI_Recv(&resource, 1, MPI_INT, MANAGER, 3, MPI_COMM_WORLD, &status);
		if (resource == -1)
		{
			MPI_Finalize();
			exit(0);
		}
		if (resource != 0)
		{ 
			resources.push(resource);
			cout << "Consumer " << i1.rank << ": i`ve got resource " << resource << endl;
			resources_to_consume--; 
		}
	}

public:
	Consumer(int in_rank, int in_resource_num) { 
		resources_to_consume = in_resource_num; 
		i1.rank = in_rank; 
		i1.whatDoYouNeed = GET_RESOURCE;
	}
	void Run() {
		while (resources_to_consume) { 
			RequestResource(); 
			RecieveResource();
		}
		i1.whatDoYouNeed = EXITC;
		MPI_Send(&i1, 3, MPI_INT, MANAGER, 0, MPI_COMM_WORLD);

	}
};

class Producer {
private:
	int resources_to_produce; 
	std::queue<int> resources;
	info i2;
	MPI_Status status;
private:
	void SendResourceToManager() { 
		i2.res = resources.front();
		int otvet;
		MPI_Send(&i2, 3, MPI_INT, MANAGER, 0, MPI_COMM_WORLD);
		
		MPI_Recv(&otvet, 1, MPI_INT, MANAGER, 1, MPI_COMM_WORLD, &status);
		if (otvet == -1)
		{
			MPI_Finalize();
			exit(0);
		}
		if (otvet == 0) 
		{
			cout << "Producer " << i2.rank << ": sending resource " << i2.res << " in buffer" << endl; 
			resources.pop(); 
		}
		else
		{
			cout << "Producer " << i2.rank << ": buffer is full. Failed to put the resource " << i2.res << " in buffer" << endl;
		}
	}
	void CreateResource() { 
		int resource = resources_to_produce + (i2.rank - 1) * 5 + 1;
		resources.push(resource);
		resources_to_produce--; 
	}
public:
	Producer(int in_rank, int num) { 
		i2.rank = in_rank; 
		resources_to_produce = num;  
		i2.whatDoYouNeed = PUT_RESOURCE;

	}
	void Run()
	{ //запуск
		while (resources_to_produce) {

			CreateResource();
			SendResourceToManager(); 
		}
		while (!resources.empty()) { 
			SendResourceToManager(); 
		}
		i2.whatDoYouNeed = EXITP;
		MPI_Send(&i2, 3, MPI_INT, MANAGER, 0, MPI_COMM_WORLD);
	}
};

class Manager { 
private:
	int total_resources; 

	int* buffer;
	int N, p_size, pro, con;
	MPI_Status status;
	info i3;

private:
	void Put(int producer_id, int resource) { 
		int otvet = 1;
		for (int i = 0; i < N; i++)
		{
			if (buffer[i] == 0) {
				otvet = 0;
				buffer[i] = resource;
				break;
			}


		}
		if (otvet == 0)
		{
			cout << "Manager: producer " << producer_id << " puts resource " << resource << endl; 
		}
		MPI_Send(&otvet, 1, MPI_INT, producer_id, 1, MPI_COMM_WORLD); 
	}
	void Get(int consumer_id) { 
		int resource = 0;
		for (int i = 0; i < N; i++)
		{
			if (buffer[i] == 0)
				resource = 0;
			else
			{
				resource = buffer[i];
				buffer[i] = 0;
				break;
			}
		}

		if (resource != 0)
		{
			cout << "Manager: consumer " << consumer_id << " gets resource " << resource << endl; 
		}
		else
		{
			if (pro == 0)
			{
				for (int i = 1; i < p_size; i++)
				{

					if (i % 2)
						MPI_Send(&STOP, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
					else
						MPI_Send(&STOP, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
				}
			}
		}
		MPI_Send(&resource, 1, MPI_INT, consumer_id, 3, MPI_COMM_WORLD);
	}

public:
	Manager(int in_total_resources, int proc_size, int p, int c) { 
		total_resources = in_total_resources;
		p_size = proc_size;
		N = total_resources;
		con = c;
		pro = p;
		buffer = new int[total_resources];
		for (int i = 0; i < total_resources; i++)
		{
			buffer[i] = 0;
		}
	}

	void Run() {
		while (true) {
			MPI_Recv(&i3, 3, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			if (i3.whatDoYouNeed == EXITP)
			{
				pro--;
			}
			if (i3.whatDoYouNeed == EXITC)
			{
				con--;
				if (con == 0)
				{
					for (int i = 1; i < p_size; i++)
					{
						if (i % 2)
							MPI_Send(&STOP, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
						else
							MPI_Send(&STOP, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
					}
					MPI_Finalize();
					exit(0);
				}
			}

			if (i3.whatDoYouNeed == PUT_RESOURCE)
				Put(i3.rank, i3.res);
			if (i3.whatDoYouNeed == GET_RESOURCE)
				Get(i3.rank);
		}
	}
};

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);
	MPI_Status status;

	int rank = -1;
	int process_num = -1;
	int pro, com;
	int size = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &process_num); 

	int a = (process_num - 1);
	pro = a / 2 + a % 2;
	com = a / 2 ;

	if (rank == 0) { 
		Manager manager(3, process_num, pro, com);
		manager.Run(); 
	}
	else { 
		if (rank % 2) { 
			std::cout << "The process with the rank of " << rank << " is Producer" << std::endl;
			Producer producer(rank, 5);
			producer.Run();
		}
		else { 
			std::cout << "The process with the rank of " << rank << " is Consumer" << std::endl;
			Consumer consumer(rank, 5); 
			consumer.Run(); 
		}
	}

	MPI_Finalize();
	return 0;
}
