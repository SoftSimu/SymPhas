
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 * 
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 * 
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 * 
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 */


#include "colormap.h"
#pragma once
#pragma warning(push, 0)   


#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkUniformGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkProgrammableFilter.h>
#include <vtkCallbackCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkSimpleElevationFilter.h>
#include <vtkPlaneSource.h>
#include <vtkNew.h>
#include <vtkColorSeries.h>
#include <vtkPolyDataMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkCommand.h>
#include <vtkAlgorithmOutput.h>
#include <vtkProgrammableFilter.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapToColors.h>
#include <vtkImageMapper3D.h>
#include <vtkImageImport.h>


#pragma warning(pop)

#include <thread>

struct filter_params
{
	filter_params(vtkImageImport* imageImport) : imageImport{ imageImport }, m{}, render_flag{ false } {}

	vtkImageImport* imageImport;
	
	std::mutex m;
	std::condition_variable cv;
	bool render_flag;
};


void timercallback(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* clientData, void* vtkNotUsed(callData))
{
	//vtkSmartPointer<vtkDoubleArray> scalars = static_cast<vtkDoubleArray*>(clientData);
	filter_params& p = *static_cast<filter_params*>(clientData);
	
	std::lock_guard<std::mutex> lock(p.m);
	if (p.render_flag)
	{
		vtkSmartPointer<vtkRenderWindowInteractor> iren = static_cast<vtkRenderWindowInteractor*>(caller);
		iren->Render();
		p.render_flag = false;
		p.cv.notify_all();

	}
}



namespace LookupTables
{

	const int red[]{ 247, 22, 22 };
	const int purple[]{ 161, 52, 235 };
	const int blue[]{ 20, 43, 222 };


	auto RedBlue(scalar_t minimum, scalar_t maximum)
	{
		// Map the scalar values in the image to colors with a lookup table:
		vtkNew<vtkLookupTable> lookupTable;
		//lookupTable->SetRange(-0.5, 0.5);
		lookupTable->SetHueRange(.05, .7);
		lookupTable->SetSaturationRange(.7, .7);
		lookupTable->SetValueRange(.6, .6);
		lookupTable->SetRampToLinear();

		lookupTable->SetNumberOfTableValues(512);
		lookupTable->SetTableRange(minimum, maximum);
		lookupTable->IndexedLookupOff();
		lookupTable->Build();

		return lookupTable;
	}

	auto OneTwoThree(scalar_t minimum, scalar_t maximum, const double one[3], const double two[2], const double three[3])
	{
		vtkNew<vtkColorTransferFunction> colorFunction;
		colorFunction->SetColorSpaceToRGB();
		colorFunction->SetColorSpaceToDiverging();
		colorFunction->AddRGBPoint(0.0, one[0], one[1], one[2]);
		colorFunction->AddRGBPoint(0.5, two[0], two[1], two[2]);
		colorFunction->AddRGBPoint(1.0, three[0], three[1], three[2]);

		vtkNew<vtkLookupTable> lookupTable;
		lookupTable->SetNumberOfTableValues(512);
		lookupTable->SetTableRange(minimum, maximum);
		lookupTable->Build();

		size_t num_colors = lookupTable->GetNumberOfTableValues();
		for (auto i = 0; i < num_colors; ++i)
		{
			double rgb[4] = { 0.0, 0.0, 0.0, 1.0 };
			colorFunction->GetColor(static_cast<double>(i) / num_colors, rgb);
			lookupTable->SetTableValue(i, rgb);
		}

		return lookupTable;
	}

	auto OneTwoThree(scalar_t minimum, scalar_t maximum, const int one[3], const int two[2], const int three[3])
	{
		double one_d[]{ one[0] / 256., one[1] / 256., one[2] / 256. };
		double two_d[]{ two[0] / 256., two[1] / 256., two[2] / 256. };
		double three_d[]{ three[0] / 256., three[1] / 256., three[2] / 256. };
		return OneTwoThree(minimum, maximum, one_d, two_d, three_d);
	}


	auto OneTwoThree(scalar_t minimum, scalar_t maximum, vtkColor3ub cone, vtkColor3ub ctwo, vtkColor3ub cthree)
	{
		int one[]{ cone.GetRed(), cone.GetGreen(), cone.GetBlue() };
		int two[]{ ctwo.GetRed(), ctwo.GetGreen(), ctwo.GetBlue() };
		int three[]{ cthree.GetRed(), cthree.GetGreen(), cthree.GetBlue() };
		return OneTwoThree(minimum, maximum, one, two, three);
	}

	auto RedPurpleBlue(scalar_t minimum, scalar_t maximum)
	{
		return OneTwoThree(minimum, maximum, red, purple, blue);
	}
	auto ColorSeriesCopy(scalar_t minimum, scalar_t maximum, vtkColorSeries::ColorSchemes cs)
	{
		vtkNew<vtkColorSeries> colorSeries;
		colorSeries->SetColorScheme(cs);


		vtkNew<vtkColorTransferFunction> colorFunction;
		colorFunction->SetColorSpaceToRGB();
		colorFunction->SetColorSpaceToDiverging();

		for (iter_type i = 0; i < colorSeries->GetNumberOfColors(); ++i)
		{
			unsigned char* color_unscaled = colorSeries->GetColor(i).GetData();
			double color[3]{};
			for (iter_type n = 0; n < 3; ++n)
			{
				color[n] = color_unscaled[n] / 256.;
			}
			colorFunction->AddRGBPoint(i / static_cast<double>(colorSeries->GetNumberOfColors() - 1), color[0], color[1], color[2]);
		}


		vtkNew<vtkLookupTable> lookupTable;
		lookupTable->SetNumberOfTableValues(512);
		lookupTable->SetTableRange(minimum, maximum);
		lookupTable->Build();

		size_t num_colors = lookupTable->GetNumberOfTableValues();
		for (auto i = 0; i < num_colors; ++i)
		{
			double rgb[4] = { 0.0, 0.0, 0.0, 1.0 };
			colorFunction->GetColor(static_cast<double>(i) / num_colors, rgb);
			lookupTable->SetTableValue(i, rgb);
		}

		return lookupTable;
	}

};


struct ColourPlotUpdaterConcrete : ColourPlotUpdater
{
	ColourPlotUpdaterConcrete(filter_params *data) : data{ data } {}
	filter_params *data;

	virtual void update()
	{
		{
			std::lock_guard<std::mutex> lock(data->m);
			data->imageImport->Modified();
			data->render_flag = true;
		}

		std::unique_lock<std::mutex> lock(data->m);
		data->cv.wait(lock, [&] { return !data->render_flag; });

	}
};


void ColourPlot2d::init(scalar_t* (&values), len_type* dims, iter_type& index, ColourPlotUpdater* (&updater))
{
	static bool one = false;
	
	if (!one)
	{
		size_t len = grid::length<2>(dims);

		scalar_t minimum = DBL_MAX, maximum = -DBL_MAX;
		for (iter_type i = 0; i < dims[0] * dims[1]; ++i)
		{
			minimum = std::min(minimum, values[i]);
			maximum = std::max(maximum, values[i]);
		}
		scalar_t range = maximum - minimum;
		scalar_t average = range / 2 + minimum;
		scalar_t colour_range = std::max(0.25, range);

		minimum = average - colour_range / 2;
		maximum = average + colour_range / 2;


		vtkNew<vtkDoubleArray> scalars;
		scalars->SetArray(values, len, 1);

		vtkNew<vtkPlaneSource> plane;
		plane->SetXResolution(dims[0] - 1);
		plane->SetYResolution(dims[1] - 1);
		plane->SetOrigin(0, 0, 0);
		plane->SetPoint1(1.0, 0.0, 0);
		plane->SetPoint2(0.0, -static_cast<double>(dims[1]) / dims[0], 0);
		plane->Update();
		plane->GetOutput()->GetPointData()->SetScalars(scalars);




		vtkNew<vtkImageData> imageData;
		imageData->SetDimensions(dims[0], dims[1], 1);
		imageData->GetPointData()->SetScalars(scalars);

		vtkNew<vtkImageImport> imageImport;
		imageImport->SetDataSpacing(.1, .1, .1);
		imageImport->SetDataOrigin(0, 0, 0);
		imageImport->SetWholeExtent(0, dims[0] - 1, 0, dims[1] - 1, 0, 0);
		imageImport->SetDataExtentToWholeExtent();
		imageImport->SetDataScalarTypeToDouble();
		imageImport->SetNumberOfScalarComponents(1);
		imageImport->SetImportVoidPointer(values);
		imageImport->Update();


		vtkNew<vtkImageMapToColors> mapColors;
		mapColors->SetLookupTable(LookupTables::RedBlue(minimum, maximum));
		mapColors->SetInputConnection(imageImport->GetOutputPort());


		vtkNew<vtkImageActor> imageActor;
		imageActor->GetMapper()->SetInputConnection(mapColors->GetOutputPort());

		//vtkNew<vtkSimpleElevationFilter> filter;
		//filter->SetInputConnection(image->GetOutputPort());
		//filter->SetVector(1, 1, 1);

		//scalarValuesToColors->SetLookupTable(lookupTable);
		//scalarValuesToColors->PassAlphaToOutputOn();
		//scalarValuesToColors->SetInputData(image);

		vtkNew<vtkPolyDataMapper> mapper;
		mapper->SetInputConnection(plane->GetOutputPort());

		mapper->SetScalarRange(minimum, maximum);
		//mapper->ScalarVisibilityOn();
		mapper->SetLookupTable(LookupTables::RedBlue(minimum, maximum));
		//plotMapper->GetProperty()->SetInterpolationTypeToNearest();

		vtkNew<vtkActor> plotActor;
		plotActor->SetMapper(mapper);

		vtkNew<vtkRenderer> renderer;
		vtkNew<vtkRenderWindow> renderWindow;
		vtkNew<vtkInteractorStyleImage> style;
		vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;

		renderWindow->SetSize(1000, 500);
		renderWindow->SetWindowName("SymPhas");


		try
		{
			renderWindow->AddRenderer(renderer);
			renderWindowInteractor->SetInteractorStyle(style);
			renderWindowInteractor->SetRenderWindow(renderWindow);

			// first initialize the render interactor
			renderWindowInteractor->Initialize();
		}
		catch (...) {}

		filter_params* p = new filter_params{ imageImport };

		// now we can add a timer based on the callback method
		vtkNew<vtkCallbackCommand> timerCallback;
		timerCallback->SetCallback(timercallback);
		timerCallback->SetClientData(p);

		renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);
		renderWindowInteractor->CreateRepeatingTimer(10);

		// Visualize
		renderer->AddActor(imageActor);
		renderer->ResetCamera();

		one = true;
		updater = new ColourPlotUpdaterConcrete{ p };

		renderWindow->Render();
		renderWindowInteractor->Start();
	}
	else
	{
		updater = new ColourPlotUpdater{};
	}
}





