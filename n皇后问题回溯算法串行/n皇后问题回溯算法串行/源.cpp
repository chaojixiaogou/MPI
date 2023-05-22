#include<iostream>
#include<cstring>
#include<Windows.h>
using namespace std;
class TimerCounter
{
public:
    TimerCounter(void);//构造函数
    ~TimerCounter(void);//析构函数
private:
    LARGE_INTEGER startCount;//记录开始时间

    LARGE_INTEGER endCount;//记录结束时间

    LARGE_INTEGER freq;//本机CPU时钟频率
public:
    double dbTime;//程序运行的时间保存在这里

public:
    void Start();//被测程序开始点处开始计时
    void Stop();//被测程序结束点处结束计时
};
TimerCounter::TimerCounter(void)
{
    QueryPerformanceFrequency(&freq);//获取主机CPU时钟频率
}
TimerCounter::~TimerCounter(void)
{
}
void TimerCounter::Start()
{
    QueryPerformanceCounter(&startCount);//开始计时
}
void TimerCounter::Stop()
{
    QueryPerformanceCounter(&endCount);//停止计时

    dbTime = ((double)endCount.QuadPart - (double)startCount.QuadPart) / (double)freq.QuadPart;//获取时间差

}

int n; // 棋盘大小
int cnt = 0; // 解的数量
int row[20]; // 存放每一行皇后的位置

// 判断在第x行放置皇后是否合法
bool check(int x) {
    for (int i = 1; i < x; i++) {
        if (row[i] == row[x] || abs(row[x] - row[i]) == abs(x - i))
            return false;
    }
    return true;
}

// 递归函数，放置第x行的皇后
void dfs(int x) {
    if (x > n) {
        cnt++;
        return;
    }
    for (int i = 1; i <= n; i++) {
        row[x] = i;
        if (check(x)) {
            dfs(x + 1);
        }
    }
}

int main() {
    cin >> n;;
    TimerCounter tc;
    tc.Start();
    dfs(1);
    tc.Stop();
    cout << cnt << endl;
    cout << tc.dbTime;
    return 0;
}
