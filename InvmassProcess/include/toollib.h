// Tools defined for Invmass process 
// Wrote by ZhengqingWang 2024.08.26 
// mail: zhengqing.wang@cern.ch
#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TProfile.h>
#include "string"
#include <TLatex.h>
#include <TLine.h>

#ifndef _TOOLLIB_H_
#define _TOOLLIB_H_
using namespace std;

class tool
{
  public:
    tool();
    ~tool();
    static TFile* FileReader(const TString &filename);
    static void SetTreeBranchAddress(TTree *tree, std::vector<Int_t*>BranchAddress_int, std::vector<Float_t*>BranchAddress_float, std::vector<Double_t*>BranchAddress_double, std::vector<TBranch*>Branches);
    static TString create_folder(const TString& path);
    static TCanvas* GetCanvas(TString title,int xpos,int ypos,int sizeX,int sizeY,Bool_t gridx,Bool_t gridy,float topMgn,float botMgn,float leftMgn,float rightMgn);
    static TPad* GetPad(TString name,float xpos1,float ypos1,float xpos2,float ypos2,float topMar,float botMar,float leftMar,float rightMar);
    static void SetMarkerTH1(TH1 *h1,TString hTitle,int markSyle,float markSize,int markColor,int lineColor);
    static void SetMarkerTH1(TGraphErrors *h1,TString hTitle,int markSyle,float markSize,int markColor,int lineColor);
    static void SetMarkerTH1(TH2 *h1,TString hTitle,int markSyle,float markSize,int markColor,int lineColor);
    static void SetTitleTH1(TH1 *h1,TString yTitle,float yTileSize,float yOffset,TString xTitle,float xTileSize,float xOffset);
    static void SetTitleTH1(TH2 *h1,TString yTitle,float yTileSize,float yOffset,TString xTitle,float xTileSize,float xOffset);
    static void SetAxisTH1(TH1 *h1,float yAxisLow,float yAxisHigh,float xAxisLow,float xAxisHigh,float yLabelsize,float xLabelsize);
    static void SetAxisTH1(TH2 *h1,float yAxisLow,float yAxisHigh,float xAxisLow,float xAxisHigh,float yLabelsize,float xLabelsize);
    static void SetAxisTH1(TH2 *h1,float yAxisLow,float yAxisHigh,float xAxisLow,float xAxisHigh,float yLabelsize,float xLabelsize,float zLabelsize,float zLabOffset);
    static void drawMyline(Float_t xstrt,Float_t ystrt,Float_t xend,Float_t yend,Int_t iStyle,Int_t iWidth,Int_t icol);
    static void drawMyText(Float_t xPos, Float_t yPos, Float_t size, TString text);
    static void drawMyTextNDC(Float_t xPos, Float_t yPos, Float_t size, TString text);//在画布的NDC(画布的归一化坐标)坐标下画
    static void drawMyTextNDC(Float_t xPos, Float_t yPos, Float_t size, TString text,Color_t col);
    static void drawBox(Float_t xstart,Float_t ystart,Float_t xstop,Float_t ystop, Int_t color, Int_t style);
};

#endif
//marker style list
//20: 实心圆点 (Full Circle)
//21: 实心方形 (Full Square)
//22: 实心三角形朝上 (Full Up-Triangle)
//23: 实心三角形朝下 (Full Down-Triangle)
//24: 实心菱形 (Full Diamond)
//25: 实心十字 (Full Cross)
//26: 实心星形 (Full Star)
//27: 空心圆点 (Open Circle)
//28: 空心方形 (Open Square)
//29: 空心三角形朝上 (Open Up-Triangle)
//30: 空心三角形朝下 (Open Down-Triangle)
//31: 空心菱形 (Open Diamond)
//32: 空心十字 (Open Cross)
//33: 空心星形 (Open Star)
//34: 星形 (Star)
//35: 星形加 (Plus)
//36: 星形 x (X Cross)
//37: 星形圆点 (Circle Star)
//=====================================================
//marker color list
//1: 黑色 (Black)
//2: 红色 (Red)
//3: 绿色 (Green)
//4: 蓝色 (Blue)
//5: 黄色 (Yellow)
//6: 粉色 (Magenta)
//7: 青色 (Cyan)
//8: 深蓝 (Dark Blue)
//9: 深绿 (Dark Green)
//10: 浅蓝 (Light Blue)
//11: 浅绿 (Light Green)
//12: 浅红 (Light Red)
//13: 灰色 (Gray)
//=====================================================
