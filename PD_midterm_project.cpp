//程式設計期中報告
//組名：白熊


#include <iostream>
#include <cmath>
using namespace std;

int main() {
    //第一大段，輸入
    int n, m, c;
    cin >> n >> m >> c;
    int nX[50] = {0};
    int nY[50] = {0};
    
    for(int i = 0; i < n; i++) {
        cin >> nX[i] >> nY[i];
    }
    
    int mX[1000] = {0};
    int mY[1000] = {0};
    
    for(int i = 0; i < m; i++) {
        cin >> mX[i] >> mY[i];
    }
    
    int mD[1000] = {0};
    
    for(int i = 0; i < m; i++) {
        cin >> mD[i];
    }
    
    int mF[1000] = {0};
    
    for(int i = 0; i < m; i++) {
        cin >> mF[i];
    }
    
    int mP[1000] = {0};
    
    for(int i = 0; i < m; i++) {
        cin >> mP[i];
    }
    
    int nH[50] = {0};
    
    for(int i = 0; i < n; i++) {
        cin >> nH[i];
    }
    
    int nK[50] = {0};
    
    for(int j = 0; j < n; j++) {
        cin >> nK[j];
    }
    
    int mDOrgn [1000] = {0};
    int nKOrgn [50] = {0};
    
    for(int i = 0; i < m; i++) {
        mDOrgn[i] = mD[i];
    }
    
    for(int j = 0; j < n; j++) {
        nKOrgn[j] = nK[j];
    }
    
    int ntom[1000][50] = {0};
    
    double nCostPer[50] = {0};
    
    for(int j = 0; j < n; j++) {
        nCostPer[j] = nH[j] / nK[j];
    }
    
    double mCostPer[1000] = {0};
    
    for(int i = 0; i < m; i++) {
        mCostPer[i] = mF[i] / mD[i];
    }
    
    //第二大段，找出每單位淨利最高的物流中心
    int jtoi[50] = {0};
    
    while(true) {
        int jMaxIndex = -1;
        int jMaxProfit = -1;
        for(int j = 0; j < n; j++) {
            if(jtoi[j] == -1) {
                continue;
            }
            int itoj[1000] = {0};
            int jProfit = 0;
            while(true) {
                int iMaxIndex = -1;
                double iMaxProfit = -1;
                for(int i = 0; i < m; i++) {
                    if(itoj[i] == -1) {
                        continue;
                    }
                    if(mD[i] == 0) {
                        itoj[i] = -1;
                        continue;
                    }
                    int d = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                    double profit = mD[i] * (mP[i] - d * c - nCostPer[j]) - mF[i];
                    profit /= mD[i];
                    if(profit > 0) {
                        if(profit > iMaxProfit) {
                            iMaxProfit = profit;
                            iMaxIndex = i;
                        }
                    } else {
                        itoj[i] = -1;
                    }
                }
                if(iMaxIndex == -1) {
                    break;
                }
                
                int i = iMaxIndex;
                int d = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                
                if(nK[j] >= mD[i]) {
                    jProfit = jProfit + mD[i] * (mP[i] - d * c) - mF[i];
                    itoj[i] = -1;
                    nK[j] -= mD[i];
                } else {
                    jProfit = jProfit + nK[j] * (mP[i] - d * c) - mF[i];
                    nK[j] = 0;
                }
                if(nK[j] == 0) {
                    break;
                }
            }
            
            if(nK[j] == nKOrgn[j]) {
                jtoi[j] = -1;
                continue;
            }
            jProfit -= nH[j];
            jProfit = jProfit / (nKOrgn[j] - nK[j]);
            if(jProfit > 0) {
                if(jProfit > jMaxProfit) {
                    jMaxProfit = jProfit;
                    jMaxIndex = j;
                }
            } else {
                jtoi[j] = -1;
            }
            nK[j] = nKOrgn[j];
        }
        
        if(jMaxIndex == -1) {
            break;
        }
        
        int j = jMaxIndex;
        jtoi[j] = -1;
        
        int itoj[1000] = {0};
        while(true) {
            int iMaxIndex = -1;
            double iMaxProfit = -1;
            for(int i = 0; i < m; i++) {
                if(itoj[i] == -1) {
                    continue;
                }
                if(mD[i] == 0) {
                    itoj[i] = -1;
                    continue;
                }
                int d = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                double profit = mD[i] * (mP[i] - d * c - nCostPer[j]) - mF[i];
                profit /= mD[i];
                if(profit > 0) {
                    if(profit > iMaxProfit) {
                        iMaxProfit = profit;
                        iMaxIndex = i;
                    }
                } else {
                    itoj[i] = -1;
                }
            }
            if(iMaxIndex == -1) {
                break;
            }
            
            int i = iMaxIndex;
            
            if(nK[j] >= mD[i]) {
                itoj[i] = -1;
                ntom[i][j] = mD[i];
                nK[j] -= mD[i];
                mD[i] = 0;
            } else {
                ntom[i][j] = nK[j];
                mD[i] -= nK[j];
                nK[j] = 0;
            }
            if(nK[j] == 0) {
                break;
            }
            
        }
        
    }
    
    //第三大段，將已選定的物流中心獲利最大化，並拿掉獲利中途變負的物流中心
    
    bool nS2[50] = {false};
    bool mS2[1000] = {false};
    
    for(int i = 0; i < m ; i++) {
        for(int j = 0; j < n; j++) {
            if(ntom[i][j] != 0) {
                nS2[j] = true;
                mS2[i] = true;
            }
        }
    }
    bool vicCtr = true;
    while(vicCtr) {
        vicCtr = false;
        //若有能使總獲利變大的零售店跟物流中心組合，就交換，或是補上
        for(int ii = 0; ii < m; ii++) {
            if(mD[ii] == 0) {
                continue;
            }
            for(int j = 0; j < n; j++) {
                if(nS2[j] == false) {
                    continue;
                }
                if(nK[j] > 0) {
                    int d = abs(mX[ii] - nX[j]) + abs(mY[ii] - nY[j]);
                    int profit = 0;
                    if(mD[ii] == mDOrgn[ii]) {
                        profit = (mP[ii] - d * c) * mD[ii] - mF[ii];
                    } else {
                        profit = (mP[ii] - d * c) * mD[ii];
                    }
                    
                    if(profit > 0) {
                        vicCtr = true;
                        if(nK[j] >= mD[ii]) {
                            ntom[ii][j] += mD[ii];
                            nK[j] -= mD[ii];
                            mD[ii] = 0;
                            break;
                        } else {
                            ntom[ii][j] += nK[j];
                            mD[ii] -= nK[j];
                            nK[j] = 0;//j-- continue
                        }
                    }
                    
                } else if (nK[j] == 0) {
                    for(int i = 0; i < m; i++) {
                        if(ii == i) {
                            continue;
                        }
                        if(ntom[i][j] == 0) {
                            continue;
                        }
                        int d = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                        int dii = abs(mX[ii] - nX[j]) + abs(mY[ii] - nY[j]);
                        if(mD[ii] == mDOrgn[ii]) {
                            if(mD[ii] >= ntom[i][j]) {
                                int profitii = (mP[ii] - dii * c) * ntom[i][j] - mF[ii];
                                int profit = (mP[i] - d * c) * ntom[i][j];
                                if(profitii > profit) {
                                    ntom[ii][j] += ntom[i][j];
                                    mD[i] += ntom[i][j];
                                    mD[ii] -= ntom[i][j];
                                    ntom[i][j] = 0;
                                    vicCtr = true;
                                }
                            } else {
                                int profitii = (mP[ii] - dii * c) * mD[ii] - mF[ii];
                                int profit = (mP[i] - d * c) * mD[ii];
                                if(profitii > profit) {
                                    ntom[ii][j] += mD[ii];
                                    mD[i] += mD[ii];
                                    ntom[i][j] -= mD[ii];
                                    mD[ii] = 0;
                                    vicCtr = true;
                                }
                            }
                        } else if (mD[ii] != mDOrgn[ii]) {
                            if(mD[ii] >= ntom[i][j]) {
                                int profitii = (mP[ii] - dii * c) * ntom[i][j];
                                int profit = (mP[i] - d * c) * ntom[i][j];
                                if(profitii > profit) {
                                    ntom[ii][j] += ntom[i][j];
                                    mD[i] += ntom[i][j];
                                    mD[ii] -= ntom[i][j];
                                    ntom[i][j] = 0;
                                    vicCtr = true;
                                }
                            } else {
                                int profitii = (mP[ii] - dii * c) * mD[ii] ;
                                int profit = (mP[i] - d * c) * mD[ii];
                                if(profitii > profit) {
                                    ntom[ii][j] += mD[ii];
                                    mD[i] += mD[ii];
                                    ntom[i][j] -= mD[ii];
                                    mD[ii] = 0;
                                    vicCtr = true;
                                }
                            }
                        }
                    }
                    
                }
            }
        }
        //物流中心跟物流中心比較，能使總獲利變大就交換
        for(int j = 0; j < n; j++){
            if(nS2[j] == false) {
                continue;
            }
            
            for(int i = 0; i < m; i++) {
                if(ntom[i][j] == 0) {
                    continue;
                }
                for(int jj = 0; jj < n; jj++) {
                    if(jj == j) {
                        continue;
                    }
                    
                    if(nS2[jj] == false) {
                        continue;
                    }
                    
                    if(nK[jj] > 0) {
                        int d = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                        int djj = abs(mX[i] - nX[jj]) + abs(mY[i] - nY[jj]);
                        if(djj < d) {
                            vicCtr = true;
                            if(ntom[i][j] >= nK[jj]) {
                                ntom[i][jj] += nK[jj];
                                ntom[i][j] -= nK[jj];
                                nK[j] += nK[jj];
                                nK[jj] = 0;
                            } else {
                                ntom[i][jj] += ntom[i][j];
                                nK[j] += ntom[i][j];
                                nK[jj] -= ntom[i][j];
                                ntom[i][j] = 0;
                                break;
                            }
                        }//else if djj < d
                        
                    }
                    
                    for(int ii = 0; ii < m; ii++) {
                        if(ntom[ii][jj] == 0) {
                            continue;
                        }
                        int d = abs(mX[ii] - nX[j]) + abs(mY[ii] - nY[j]);
                        int djj = abs(mX[i] - nX[jj]) + abs(mY[i] - nY[jj]);
                        int proftiPer = mP[ii] - d * c;
                        int profitPerjj = mP[i] - djj * c;
                        
                        int dOrg = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                        int dOrgjj = abs(mX[ii] - nX[jj]) + abs(mY[ii] - nY[jj]);
                        int profitperOrg = mP[i] - dOrg * c;
                        int profitperOrgjj = mP[ii] - dOrgjj * c;
                        
                        if(proftiPer + profitPerjj - profitperOrg - profitperOrgjj > 0) {
                            vicCtr = true;
                            if(ntom[i][j] >= ntom[ii][jj]) {
                                ntom[ii][j] += ntom[ii][jj];
                                ntom[i][jj] += ntom[ii][jj];
                                ntom[i][j] -= ntom[ii][jj];
                                ntom[ii][jj] = 0;
                            } else if(ntom[ii][jj] > ntom[i][j]) {
                                ntom[ii][j] += ntom[i][j];
                                ntom[i][jj] += ntom[i][j];
                                ntom[ii][jj] -= ntom[i][j];
                                ntom[i][j] = 0;
                            }
                        }
                    }
                }
            }
        }
        //拿掉獲利變負的物流中心
        for(int j = 0; j < n; j++) {
            if(nS2[j] == false) {
                continue;
            }
            int jTotalProfit = 0;
            for(int i = 0; i < m; i++) {
                if(ntom[i][j] != 0) {
                    int d = abs(mX[i] - nX[j]) + abs(mY[i] - nY[j]);
                    jTotalProfit = jTotalProfit + ntom[i][j] * (mP[i] - d * c) - mF[i];
                }
            }
            jTotalProfit -= nH[j];
            if(jTotalProfit < 0) {
                for(int i = 0; i < m; i++) {
                    if(ntom[i][j] != 0) {
                        mD[i] += ntom[i][j];
                        nK[j] += ntom[i][j];
                        ntom[i][j] = 0;
                    }
                }
                nS2[j] = false;
                vicCtr = true;
            }
        }
        
        
    }
    
    //大四大段，debug跟輸出
    bool mSOut[1000] = {false};
    bool nSOut[50] = {false};
    
    for(int j = 0; j < n; j++) {
        if(nK[j] < 0) {
            return 0;
        }
        int jTest = 0;
        for(int i = 0; i < m; i++) {
            if(ntom[i][j] < 0 || mD[i] < 0) {
                return 0;
            } else {
                jTest += ntom[i][j];
            }
        }
        jTest += nK[j];
        if(jTest != nKOrgn[j]) {
            return 0;
        }
    }
    
    for(int i = 0; i < m; i++) {
        int iTest = 0;
        for(int j = 0; j < n; j++) {
            if(ntom[i][j] != 0) {
                iTest += ntom[i][j];
            }
        }
        iTest += mD[i];
        if(iTest != mDOrgn[i]) {
            return 0;
        }
    }
    
    for(int i = 0; i < m ; i++) {
        for(int j = 0; j < n; j++) {
            if(ntom[i][j] != 0) { //改!= 0試試
                mSOut[i] = true;
                nSOut[j] = true;
            }
        }
    }
    
    int nCnt = 0, mCnt = 0;
    
    for(int j = 0; j < n; j++) {
        if(nSOut[j]) {
            nCnt++;
        }
    }
    
    for(int i = 0; i < m; i++) {
        if(mSOut[i]) {
            mCnt++;
        }
    }
    
    cout << nCnt;
    
    for(int j = 0; j < n; j++) {
        if(nSOut[j]) {
            cout << " " << j + 1;
        }
    }
    
    cout << endl << mCnt;
    
    for(int i = 0; i < m; i++) {
        if(mSOut[i]) {
            cout << " " << i + 1;
        }
    }
    
    cout << endl;
    
    int a = 0, b = 0;
    
    for(int i = 0; i < m; i++) {
        a = 0;
        cout << (b==0?"":"\n");
        
        for(int j = 0; j < n; j++) {
            cout << (a==0?"":" ") << ntom[i][j];
            a = 1;
        }
        b = 1;
    }
    
    return 0;
}
