/*******************************************************************

VBTest.cpp					developed by naka_t	2011.01.31

			変分ベイズGMMクラスの使用例

Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
*******************************************************************/
#include "stdafx.h"
#include "VariationalBayes.h"
#include "utility.h"

int _tmain(int argc, _TCHAR* argv[])
{
	int dataNum = 0;
	int dim = 0;
	double **data = NULL;
	CVariationalBayes vb;

	// データの読み込み
	data = LoadMatrix<double>( dim , dataNum , "data.txt" );

	vb.SetData( data , dataNum );	// データを入力
	vb.Learn( "result" );			// 学習開始

	printf("最適混合数 %d\n" , vb.GetBestModelNum() );

	Free( data );
	return 0;
}

