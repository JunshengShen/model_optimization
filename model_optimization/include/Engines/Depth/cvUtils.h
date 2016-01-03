#ifndef CV_UTILS_H_
#define CV_UTILS_H_

#include  <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/ml/ml.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace SU{
	class suCvUtils{
	public:

		suCvUtils();

		//utilities functions
		static void writeToTxt(char*pfilename, cv::Mat1s m){
			std::ofstream f(pfilename,std::ios::out);
			for (int i=0;i<m.rows; i++)
			{
				for(int j=0; j<m.cols; j++)
				{
					f << m.at<short>(i,j) << ", ";
				}
				f << "\n";
			}
			f.close();
		}

/*
 * \brief readDepthFromTxt
 * \parameters
 *@strFile file name
 *@outMat cv::Mat dMat(rows,cols, CV_16UC1);
 */
		static void readDepthFromTxt( std::string strFile, int rows, int cols, cv::Mat &outMat)
		{
			std::ifstream infile(strFile.c_str());
			std::vector <std::vector <std::string> > data;
			int h=rows,w=cols;
			while (infile)
			{
				string s;
				if (!getline( infile, s )) break;

				istringstream ss( s );
				vector <string> record;

				while (ss)
				{
					string s;
					if (!getline( ss, s, ',' )) break;
					record.push_back( s );
				}

				data.push_back( record );
			}
			if (!infile.eof())
			{
				cerr << "file in error!\n";
			}
			int nDepthWidth = w;
			for (int y = 0; y < h; y++)
			{
				unsigned short* pDepth   = (unsigned short*)outMat.data + y*nDepthWidth;

				for (int x = 0; x < nDepthWidth; x++)
				{
					*(pDepth+x) = (unsigned short) atoi(data[y][x].c_str() );
				}
			}
		}

		/*
 * \brief save depth image to txt format file.c_str 
 * \parameters
 *@strFile file name
 *@outMat cv::Mat dMat(rows,cols, CV_16UC1);
 */
static void saveDepthToTxt( std::string strFile, cv::Mat &outMat)
{
	std::ofstream outfile(strFile.c_str());
	std::vector <std::vector <std::string> > data;
	int h=outMat.rows;
	int w=outMat.cols;

	if (outfile)
	{
		for (int y = 0; y < h; y++)
		{
			unsigned short* pDepth   = (unsigned short*)outMat.data + y * w;

			for (int x = 0; x < w - 1; x++)
			{
				outfile << *(pDepth+x) << "," ;
			}
			outfile << *(pDepth + w -1) << std::endl;
		}
	}
	
	outfile.close();	
}

/** 
 *bDrawPoint = 1: draw vertex
 */
static void drawContour(cv::Mat& m, std::vector<cv::Point> contour, cv::Scalar &color, int thickness, int bDrawPoint=0 )
{
	int n = (int)contour.size() - 1;
	for (int i=0; i<n; i++ )
	{
		cv::line(m, contour[i], contour[i+1], color, thickness);
	}
	cv::line(m, contour[n], contour[0], color, thickness);

	if (bDrawPoint == 1)
	{
		for (int i=0; i<(int)contour.size(); i++ )
		{
			cv::circle(m, contour[i], 4, color);
		}
	}
		
}
/*
 *\brief read depth from 16U png file
 * \parameters
 * @strFile file name
 * @return  cv::Mat dMat(rows,cols, CV_16UC1);
 */
		static cv::Mat readDepthFromPng( std::string strFile)
		{
			return 	cv::imread(strFile, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR );
		}


	};
	/*Draw line to a image with the line normal vector 
		*(@lineParameter[0], @lineParameters[1]) 
		*and a point on line 
		*(@lineParameter[2], @lineParameters[3]) 
		*equation:  ax + bx = c
		*           a = nx = @lineParameter[0]
		*           b = ny = @lineParameter[2]
		*/	
	void drawLine(cv::Mat mat, int x1, int x2,std::vector<double> lineParameters, int colorIdx = 1, int lineWidth=1);
	void drawLineScale(cv::Mat mat, int lenFilter,int x1, int x2,std::vector<double> lineParameters, int colorIdx, int lineWidth);
	

}

#endif