#pragma once
/*******************************************************************

utility.h		developed by naka_t	2011.01.31

	�t�@�C�����o�͌n�֐�

�@��VariationalBayes.h�ɕK�v�Ȋ֐����ڐA			naka_t	2011.01.31

Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
*******************************************************************/
#include <vector>
#include <windows.h>

#pragma warning(disable: 4996)
#pragma warning(disable: 4018)

#ifndef Message
#define Message( mess ) ::MessageBox( 0 , mess , "Message" , 0 )
#endif


// �t�@�C�����o�͊֐�
static int fprintfVar( int val , FILE *fp ){ return fprintf( fp , "%d" , val ); }
static int fprintfVar( double val , FILE *fp ){ return fprintf( fp , "%lf" , val ); }
static int fprintfVar( std::string &val , FILE *fp ){ return fprintf( fp , "%s" , val.c_str() ); }
static int fprintfVar( const char *&val , FILE *fp ){ return fprintf( fp , "%s" , val ); }
static int fscanfVar( int *var , FILE *fp ){ return fscanf( fp, "%d" , var ); }
static int fscanfVar( double *var , FILE *fp ){ return fscanf( fp, "%lf" , var ); }
static int fscanfVar( float *var , FILE *fp ){ return fscanf( fp, "%f" , var ); }
static int fscanfVar( bool *var , FILE *fp ){ return fscanf( fp, "%d" , var ); }
static int fscanfVar( std::string *var , FILE *fp ){
	char tmpStr[1024] = {0};
	int ret = fscanf( fp , "%s" , tmpStr );
	*var = tmpStr;
	return ret;
}

// std::string�ł�fgets
// �߂�l�F�ǂݍ��񂾕�����
//		   �t�@�C���̏I�[�̍ۂ́Astring::npos��Ԃ��܂�
static size_t fgets( std::string &buff , FILE *fp ){
	char c;
	buff.clear();
	while ( ( c = fgetc( fp ) ) != EOF && c != '\n' )
		buff.push_back( c );

	if( c==EOF && buff.size()==0 ) return std::string::npos;

	return buff.size();
}

//���O���ŕϐ���ۑ�
template <typename T> bool SaveVariable( T value , const char *valname , const char *filename , const char *mode="a" )
{
	FILE *fp;
	std::vector<std::string> vecLines;
	std::string line;
	bool bFind = false;

	if( *mode=='a')	// �ǋL���[�h�Ȃ�Ό��݂̓��e���ꎞ�ۑ�
	{
		fp = fopen( filename , "r" );
		if( fp )	// �t�@�C��������Γǂݍ���
		{
			while( fgets( line , fp ) != std::string::npos ) { vecLines.push_back(line); }
			fclose(fp);
		}
	}


	// �������݊J�n
	fp = fopen( filename , "w" );

	if( !fp )
	{
		char mess[1024];
		sprintf( mess , "SaveVariable : %s���I�[�v���ł��܂���" , filename );
		Message( mess );
		return false;
	}

	for(int i=0 ; i<vecLines.size() ; i++ )
	{
		char tmpName[256];
		sscanf( vecLines[i].c_str() , "%s" , tmpName );
		if( strcmp(tmpName,valname)==0 )
		{
			// �ϐ�������������A�㏑��
			fprintf( fp , "%s		" , valname );
			fprintfVar( value , fp );
			fprintf( fp , "\n" );
			bFind = true;
		}
		else
		{
			// ������Ȃ�������A���̂܂ܕۑ�
			fprintf( fp , "%s\n" , vecLines[i].c_str() );
		}
	}

	if( !bFind )
	{
		fprintf( fp , "%s		" , valname );
		fprintfVar( value , fp );
		fprintf( fp , "\n" );
	}

	fclose(fp);
	return true;
}


//���O���ŕϐ���ۑ�
template <typename T> bool SaveVariableArray( T value , int num , const char *valname , const char *filename , const char *mode="a" )
{
	FILE *fp;
	std::vector<std::string> vecLines;
	std::string line;
	bool bFind = false;

	if( *mode=='a')	// �ǋL���[�h�Ȃ�Ό��݂̓��e���ꎞ�ۑ�
	{
		fp = fopen( filename , "r" );
		if( fp )	// �t�@�C��������Γǂݍ���
		{
			while( fgets( line , fp ) != std::string::npos ) { vecLines.push_back(line); }
			fclose(fp);
		}
	}


	// �������݊J�n
	fp = fopen( filename , "w" );

	if( !fp )
	{
		char mess[1024];
		sprintf( mess , "SaveVariable : %s���I�[�v���ł��܂���" , filename );
		Message( mess );
		return false;
	}

	for(int i=0 ; i<vecLines.size() ; i++ )
	{
		char tmpName[256];
		sscanf( vecLines[i].c_str() , "%s" , tmpName );
		if( strcmp(tmpName,valname)==0 )
		{
			// �ϐ�������������A�㏑��
			fprintf( fp , "%s	" , valname );
			for(int j=0 ; j<num ; j++ )
			{
				fprintfVar( value[j] , fp );
				fprintf( fp , "	" );
			}
			fprintf( fp , "\n" );
			bFind = true;
		}
		else
		{
			// ������Ȃ�������A���̂܂ܕۑ�
			fprintf( fp , "%s\n" , vecLines[i].c_str() );
		}
	}

	if( !bFind )
	{
		fprintf( fp , "%s	" , valname );
		for(int j=0 ; j<num ; j++ )
		{
			fprintfVar( value[j] , fp );
			fprintf( fp , "	" );
		}
		fprintf( fp , "\n" );
	}

	fclose(fp);
	return true;
}

//���O���ŕϐ���ۑ�
template <typename T> bool SaveVariableMatrix( T value , int ysize , int xsize , const char *valname , const char *filename , const char *mode="a" )
{
	FILE *fp;
	std::vector<std::string> vecLines;
	std::string line;
	bool bFind = false;

	if( *mode=='a')	// �ǋL���[�h�Ȃ�Ό��݂̓��e���ꎞ�ۑ�
	{
		fp = fopen( filename , "r" );
		if( fp )	// �t�@�C��������Γǂݍ���
		{
			while( fgets( line , fp ) != std::string::npos ) { vecLines.push_back(line); }
			fclose(fp);
		}
	}


	// �������݊J�n
	fp = fopen( filename , "w" );

	if( !fp )
	{
		char mess[1024];
		sprintf( mess , "SaveVariable : %s���I�[�v���ł��܂���" , filename );
		Message( mess );
		return false;
	}

	for(int i=0 ; i<vecLines.size() ; i++ )
	{
		char tmpName[256];
		sscanf( vecLines[i].c_str() , "%s" , tmpName );
		if( strcmp(tmpName,valname)==0 )
		{
			fprintf( fp , "%s			%d	%d\n" , valname , ysize , xsize );

			// �ϐ�������������A�㏑��
			for(int y=0 ; y<ysize ; y++ )
			{
				fprintf( fp , "*	" );
				for(int x=0 ; x<xsize ; x++ )
				{
					fprintfVar( value[y][x] , fp );
					fprintf( fp , "	" );
				}
				fprintf( fp , "\n" );
			}

			// "*"���s��̑�����\���̂œǂݔ�΂�
			for( i = i+1 ; i<vecLines.size() ; i++ )
			{
				sscanf( vecLines[i].c_str() , "%s" , tmpName );
				if( strcmp(tmpName,"*")!=0 ) 
				{
					i--;
					break;
				}
			}
			bFind = true;
		}
		else
		{
			// ������Ȃ�������A���̂܂ܕۑ�
			fprintf( fp , "%s\n" , vecLines[i].c_str() );
		}
	}

	if( !bFind )
	{
		fprintf( fp , "%s			%d	%d\n" , valname , ysize , xsize );
		for(int y=0 ; y<ysize ; y++ )
		{
			fprintf( fp , "*	" );
			for(int x=0 ; x<xsize ; x++ )
			{
				fprintfVar( value[y][x] , fp );
				fprintf( fp , "	" );
			}
			fprintf( fp , "\n" );
		}

	}

	fclose(fp);
	return true;
}


//�s�����J�E���g
static int CountLine( FILE *fp )
{
	int num = 0;
	int c;

	while(1)
	{
		c = fgetc(fp);

		if( c == '\0' || c == EOF ) break;
		else if( c == '\n' ) num++;
	}

	//printf( "lines %d\n" , num );

	rewind(fp);

	return num;
}

//�s���̃^�u�̐��𐔂���
static int CountCols( FILE *fp )
{
	int num = 0;
	int c , old;

	while(1)
	{
		c = fgetc(fp);

		if( c == '\n' || c == EOF ) break;
		else if( c == '\t' ) num++;

		old = c;
	}

	//�Ō�̗�Ƀ^�u�������Ă��Ȃ��\��������̂�
	if( old != '\t' && c == '\n' ) num++;

	rewind(fp);

	return num;
}

// �������m��
template<typename T> T** AllocMatrix( int ysiz , int xsiz)
{
	//�z��̎��Ԃ�1�����Ŋm�ۂ��Ă��邱�Ƃɒ��ӁI

	T **mptr = NULL;
	T *data = NULL;	//�z��̎���

	try
	{
		mptr = new T*[ysiz];
		data = new T[ysiz*xsiz];

		for (int i=0;i<ysiz;i++)
		{
			mptr[i] = data + i*xsiz;
		}

		return mptr;
	}
	catch( std::bad_alloc &ba )
	{
		std::string mess = "AllocMatrix : �������̊m�ۂɎ��s�I\r\n";
		mess += ba.what();
		Message(mess.c_str());
		return NULL;
	}
}

// ���������
template <typename T> void Free(T **mptr)
{
	if( mptr )
	{
		delete [] mptr[0];
		delete [] mptr;
	}
}

template <typename T> bool LoadMatrix( T **array , int xsize, int ysize , FILE *fp )
{
	T data;

	for(int m=0 ; m<ysize ; m++)
	{
		for(int n=0 ; n<xsize ; n++)
		{
			if( fscanfVar( &data , fp ) != 1 )
			{
				return false;
			}
			fscanf( fp , "	" );

			array[m][n]=data;
		}
		
		fscanf( fp , "\n");
	}
	return true;
}


template <typename T> T **LoadMatrix(int &xsize, int &ysize , const char *filename )
{
	FILE *fp = NULL;
	T **array = NULL;
	char mess[1024];

	fp = fopen( filename , "r");

	if( fp==NULL )
	{
		sprintf( mess, "%s ������܂���D" , filename );
		Message( mess );

	}


	//�s���擾
	ysize = CountLine(fp);
	xsize = CountCols(fp);

	array = (T **)AllocMatrix<T>( ysize , xsize );

	if( !LoadMatrix( array , xsize, ysize , fp ) )
	{
		if( !fp ) fclose( fp );
		if( !array ) Free( array );

		sprintf( mess, "%s �̏������s���ł��D" , filename );
		Message( mess );

		return NULL;
	}

	fclose(fp);

	return array;
}
