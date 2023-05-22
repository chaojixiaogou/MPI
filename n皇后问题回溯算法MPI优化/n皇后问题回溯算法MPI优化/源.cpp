#include <iostream>
#include <mpi.h>
#include<Windows.h>
using namespace std;
class TimerCounter
{
public:
    TimerCounter(void);//���캯��
    ~TimerCounter(void);//��������
private:
    LARGE_INTEGER startCount;//��¼��ʼʱ��

    LARGE_INTEGER endCount;//��¼����ʱ��

    LARGE_INTEGER freq;//����CPUʱ��Ƶ��
public:
    double dbTime;//�������е�ʱ�䱣��������

public:
    void Start();//�������ʼ�㴦��ʼ��ʱ
    void Stop();//�����������㴦������ʱ
};
TimerCounter::TimerCounter(void)
{
    QueryPerformanceFrequency(&freq);//��ȡ����CPUʱ��Ƶ��
}
TimerCounter::~TimerCounter(void)
{
}
void TimerCounter::Start()
{
    QueryPerformanceCounter(&startCount);//��ʼ��ʱ
}
void TimerCounter::Stop()
{
    QueryPerformanceCounter(&endCount);//ֹͣ��ʱ

    dbTime = ((double)endCount.QuadPart - (double)startCount.QuadPart) / (double)freq.QuadPart;//��ȡʱ���

}

int n; // ���̴�С
int cnt = 0; // �������
int row[20]; // ���ÿһ�лʺ��λ��

// �ж��ڵ�x�з��ûʺ��Ƿ�Ϸ�
bool check(int x) {
    for (int i = 1; i < x; i++) {
        if (row[i] == row[x] || abs(row[x] - row[i]) == abs(x - i))
            return false;
    }
    return true;
}

// �ݹ麯�������õ�x�еĻʺ�
void dfs(int x, int& partialCount) {
    if (x > n) {
        partialCount++;
        return;
    }
    for (int i = 1; i <= n; i++) {
        row[x] = i;
        if (check(x)) {
            dfs(x + 1, partialCount);
        }
    }
}

int main(int argc, char** argv) {
    TimerCounter tc;

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        cin >> n;
        tc.Start();
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int partialCount = 0;
    int start = rank * (n / size) + 1;
    int end = (rank + 1) * (n / size);

    if (rank == size - 1) {
        end = n;
    }

    for (int i = start; i <= end; i++) {
        row[1] = i;
        dfs(2, partialCount);
    }

    int total;
    MPI_Reduce(&partialCount, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << total << endl;
        tc.Stop();
        cout << tc.dbTime;
    }

    MPI_Finalize();


    return 0;
}
