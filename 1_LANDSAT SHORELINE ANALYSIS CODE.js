/**
 * LANDSAT SHORELINE ANALYSIS CODE - ENGLISH VERSION
 * 
 * @description
 * This script processes Landsat satellite images using Google Earth Engine 
 * to analyze surface changes in the "Sector 5" region from 1984 to 2024.
 * 
 * The script performs the following operations:
 * - Centers the map on the region of interest (ROI)
 * - Loads Landsat 4, 5, 7, 8, and 9 image collections
 * - Defines cloud mask functions for different Landsat sensors
 * - Filters and composes images applying cloud masks
 * - Calculates water indices (NDWI, MNDWI) and applies Otsu thresholding
 * - Uses Canny edge detection for coastline extraction
 * - Provides interactive UI for visualization and data export
 * - Exports results in raster and vector formats
 * 
 * @authors DIAS, K. S. & MENEZES, R. A. A.
 * @version 4.1
 * @year 2025
 */

// ============================================================================
// SECTION 1: USER INTERFACE SETUP
// ============================================================================

// Create main control panel
var panel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { width: '400px', padding: '10px' }
});

// Header and version information
var header = ui.Label({
  value: 'LANDSAT SHORELINE ANALYSIS INTERFACE',
  style: { fontSize: '18px', fontWeight: 'bold', textAlign: 'center', padding: '10px' }
});

var versionLabel = ui.Label({
  value: 'Code version: 4.1',
  style: { fontSize: '14px', textAlign: 'center', color: 'gray', padding: '5px' }
});

// Developer credits
var label = ui.Label('Developed by:');
var linkDias = ui.Label('DIAS, K. S.', { color: 'blue' }, 
                        'http://lattes.cnpq.br/5961292748412062');
var linkMenezes = ui.Label('MENEZES, R. A. A.', { color: 'blue' }, 
                        'http://lattes.cnpq.br/4681123810496065');

var developers = ui.Panel([
  label, ui.Label(' '), linkDias, ui.Label(' & '), linkMenezes
], ui.Panel.Layout.flow('horizontal'));

// Add header components to panel
panel.add(header);
panel.add(versionLabel);
panel.add(developers);

// ============================================================================
// SECTION 2: MAP INITIALIZATION AND DATA LOADING
// ============================================================================

// Center map on region of interest
Map.centerObject(roi, 8);

// Load Landsat image collections
var l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2"),
    l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2"),
    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"),
    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");

// Add Sector 5 boundary to map
var empty = ee.Image().byte();
var contour = empty.paint({
  featureCollection: table,
  color: 1,
  width: 2
});
Map.addLayer(contour, { palette: 'red' }, 'Sector 5', true);

// ============================================================================
// SECTION 3: IMAGE PROCESSING FUNCTIONS
// ============================================================================

// Cloud masking functions for different Landsat sensors
function cloudMaskTm(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 4).eq(0) // Cloud Shadow
    .and(qa.bitwiseAnd(1 << 3).eq(0));   // Clouds
  return image.updateMask(mask)
    .select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
            ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
}

function cloudMaskOli(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0) // Dilated Cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))   // Cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))   // Cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0));  // Cloud Shadow
  return image.updateMask(mask)
    .select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
            ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
}

// Filter and prepare image collections
function filterAndPrepare(col, roi, date, sensorType) {
  var maskFunction = sensorType === 'OLI' ? cloudMaskOli : cloudMaskTm;
  return col.filterBounds(roi)
            .filterDate(date[0], date[1])
            .map(maskFunction);
}

// Create Landsat composite for given year and date range
function landsatComposite(roi, date, year) {
  var col;
  if (year < 2013) {
    col = filterAndPrepare(l4, roi, date, 'TM')
      .merge(filterAndPrepare(l5, roi, date, 'TM'))
      .merge(filterAndPrepare(l7, roi, date, 'TM'));
  } else {
    col = filterAndPrepare(l8, roi, date, 'OLI')
      .merge(filterAndPrepare(l9, roi, date, 'OLI'));
  }
  var image = col.median().clip(roi);
  return image.multiply(0.0000275).add(-0.2);
}

// ============================================================================
// SECTION 4: OTSU THRESHOLDING ALGORITHM
// ============================================================================

// Calculate Otsu threshold for automatic water/land separation
function getOtsuThreshold(image, region) {
  var bandName = image.bandNames().get(0);
  var histogram = image.reduceRegion({
    reducer: ee.Reducer.histogram({ maxBuckets: 256 }),
    geometry: region,
    scale: 30,
    bestEffort: true
  });

  var histData = ee.Dictionary(histogram.get(bandName));
  print('Extracted histogram:', histData);

  var otsu = ee.Algorithms.If(
    histData.contains('histogram'), 
    ee.Number(otsuThreshold(histData)), 
    ee.Number(-999)
  );
  return otsu;
}

// Otsu algorithm implementation
function otsuThreshold(histogram) {
  var counts = ee.Array(histogram.get('histogram'));
  var means = ee.Array(histogram.get('bucketMeans'));

  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var size = counts.length().get([0]);

  // Initialize variables
  var wB = ee.Number(0);
  var sumB = ee.Number(0);
  var maxVar = ee.Number(0);
  var threshold = ee.Number(0);

  var iterate = ee.List.sequence(0, size.subtract(1)).iterate(function(i, prev) {
    i = ee.Number(i);
    prev = ee.Dictionary(prev);

    var c = counts.get([i]);
    var m = means.get([i]);

    var wB_new = ee.Number(prev.get('wB')).add(c);
    var sumB_new = ee.Number(prev.get('sumB')).add(m.multiply(c));
    var wF = ee.Number(total).subtract(wB_new);

    var mB = sumB_new.divide(wB_new);
    var mF = ee.Number(sum).subtract(sumB_new).divide(wF);
    var between = wB_new.multiply(wF).multiply(mB.subtract(mF).pow(2));

    var maxVar = ee.Number(prev.get('maxVar'));
    var threshold = ee.Number(prev.get('threshold'));

    var newThreshold = ee.Algorithms.If(between.gt(maxVar), m, threshold);
    var newMaxVar = ee.Algorithms.If(between.gt(maxVar), between, maxVar);

    return ee.Dictionary({
      wB: wB_new,
      sumB: sumB_new,
      maxVar: newMaxVar,
      threshold: newThreshold
    });
  }, {
    wB: wB,
    sumB: sumB,
    maxVar: maxVar,
    threshold: threshold
  });

  return ee.Number(ee.Dictionary(iterate).get('threshold'));
}

// ============================================================================
// SECTION 5: VISUALIZATION AND ANALYSIS FUNCTIONS
// ============================================================================

// Visualization and analysis helper functions
function plotOtsuSeparationChart(histogram, thresholdValue, year) {
  var counts = ee.Array(histogram.get('histogram'));
  var bins = ee.Array(histogram.get('bucketMeans'));
  var length = ee.Number(bins.length().get([0]));

  // Initialize class separation lists
  var aBins = ee.List([]), aFreqs = ee.List([]);
  var bBins = ee.List([]), bFreqs = ee.List([]);

  var separated = ee.List.sequence(0, length.subtract(1)).iterate(function(i, acc) {
    i = ee.Number(i);
    acc = ee.Dictionary(acc);

    var bin = ee.Number(bins.get([i]));
    var count = ee.Number(counts.get([i]));

    var aBins = ee.List(acc.get('aBins'));
    var aFreqs = ee.List(acc.get('aFreqs'));
    var bBins = ee.List(acc.get('bBins'));
    var bFreqs = ee.List(acc.get('bFreqs'));

    return ee.Algorithms.If(bin.lte(thresholdValue),
      ee.Dictionary({ 
        aBins: aBins.add(bin), aFreqs: aFreqs.add(count), 
        bBins: bBins, bFreqs: bFreqs 
      }),
      ee.Dictionary({ 
        aBins: aBins, aFreqs: aFreqs, 
        bBins: bBins.add(bin), bFreqs: bFreqs.add(count) 
      })
    );
  }, ee.Dictionary({ aBins: aBins, aFreqs: aFreqs, bBins: bBins, bFreqs: bFreqs }));

  separated = ee.Dictionary(separated);

  // Create scatter charts for each class
  var chartA = ui.Chart.array.values({
    array: ee.Array([separated.get('aFreqs')]),
    axis: 1,
    xLabels: ee.List(separated.get('aBins'))
  }).setChartType('ScatterChart').setOptions({
    title: null,
    colors: ['#000000'],
    pointSize: 5,
    hAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
    vAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
    legend: { position: 'none' },
    backgroundColor: { fill: 'transparent' }
  });

  var chartB = ui.Chart.array.values({
    array: ee.Array([separated.get('bFreqs')]),
    axis: 1,
    xLabels: ee.List(separated.get('bBins'))
  }).setChartType('ScatterChart').setOptions({
    title: null,
    colors: ['#adadad'],
    pointSize: 5,
    hAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
    vAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
    legend: { position: 'none' },
    backgroundColor: { fill: 'transparent' }
  });

  print(chartA);
  print(chartB);
}

// Plot class probabilities
function plotOtsuProbabilities(probA, probB, year) {
  var labels = ['Class A', 'Class B'];
  var values = [probA, probB];

  print(ui.Chart.array.values({
    array: ee.Array(values.map(function(v) {
      return [v];
    })),
    axis: 0,
    xLabels: labels
  }).setChartType('ScatterChart')
    .setOptions({
      title: null,
      pointSize: 6,
      series: {
        0: {color: 'black'},
        1: {color: 'gray'}
      },
      hAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
      vAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
      legend: { position: 'none' },
      backgroundColor: { fill: 'transparent' }
    }));
}

// ============================================================================
// SECTION 6: USER INTERFACE CONTROLS
// ============================================================================

// Year selection dropdown
var yearList = ee.List.sequence(1984, 2024).getInfo();
var yearDropdown = ui.Select({
  items: yearList.map(function(y) { return y.toString(); }),
  placeholder: 'Select Year'
});

// Layer visibility checkboxes
var checkboxRGB = ui.Checkbox({ label: 'Landsat Composite (RGB)', value: true });
var checkboxNDWI = ui.Checkbox({ label: 'NDWI', value: true });
var checkboxMNDWI = ui.Checkbox({ label: 'MNDWI', value: true });
var checkboxWater = ui.Checkbox({ label: 'Water', value: true });
var checkboxLand = ui.Checkbox({ label: 'Land', value: true });
var checkboxCoastline = ui.Checkbox({ label: 'Coastline', value: true });
var checkboxOtsu = ui.Checkbox({ label: 'Otsu Threshold', value: true });

// ============================================================================
// SECTION 7: MAP UPDATE AND VISUALIZATION LOGIC
// ============================================================================

// Main map update function
function updateMap() {
  var selectedYear = yearDropdown.getValue();
  if (!selectedYear) return;

  var year = parseInt(selectedYear, 10);
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year + 4, 12, 31);
  var dateRange = [startDate, endDate];
  var image = landsatComposite(roi, dateRange, year);

  Map.layers().reset();

  // RGB Visualization
  if (checkboxRGB.getValue()) {
    var visParams = { min: 0, max: 0.3 };
    if (year < 2013) {
      Map.addLayer(image.select(['B3', 'B2', 'B1']), visParams, 'Landsat_' + year); 
    } else {
      Map.addLayer(image.select(['B4', 'B3', 'B2']), visParams, 'Landsat_' + year);
    }
  }

  // Calculate water indices
  var ndwi = image.expression('(B3 - B5) / (B3 + B5)', {
    B3: image.select('B3'),
    B5: image.select('B5')
  }).rename('NDWI');

  var mndwi = image.expression('(B3 - B5) / (B3 + B5)', {
    B3: image.select('B3'),
    B5: image.select('B5')
  }).rename('MNDWI');

  // Add layers based on checkbox selections
  if (checkboxNDWI.getValue()) {
    Map.addLayer(ndwi, { min: -1, max: 1, palette: ['brown', 'white', 'blue'] }, 'NDWI_' + year);
  }

  if (checkboxMNDWI.getValue()) {
    Map.addLayer(mndwi, { min: -1, max: 1, palette: ['blue', 'white', 'green'] }, 'MNDWI_' + year);
  }

  if (checkboxWater.getValue()) {
    var water = mndwi.gte(0.1).selfMask();
    Map.addLayer(water, { palette: 'white' }, 'Water_' + year);
  }

  if (checkboxLand.getValue()) {
    var land = mndwi.lt(0.1).selfMask();
    Map.addLayer(land, { palette: 'black' }, 'Land_' + year);
  }

  if (checkboxCoastline.getValue()) {
    var coastline = ee.Algorithms.CannyEdgeDetector({
      image: mndwi, threshold: 0.7, sigma: 1
    }).selfMask();
    Map.addLayer(coastline, { palette: 'red' }, 'Coastline_' + year);
  }

  // Otsu threshold processing
  if (checkboxOtsu.getValue()) {
    var otsuThresholdValue = getOtsuThreshold(mndwi, roi);
    otsuThresholdValue.evaluate(function(thresh) {
      print('🔹 Otsu threshold for year ' + year + ':', thresh);
      if (thresh !== null && typeof thresh === 'number' && thresh !== -999) {

        // Remove previous Otsu layer if exists
        Map.layers().forEach(function(layer) {
          if (layer.getName() === 'Binary Water (Otsu) ' + year) {
            Map.layers().remove(layer);
          }
        });

        var binaryWater = mndwi.gte(thresh);

        // Statistical calculations
        var binCollection = ee.ImageCollection([binaryWater]);
        var sum = binCollection.reduce(ee.Reducer.sum());
        var count = binCollection.reduce(ee.Reducer.count());
        var fp = sum.divide(count);

        // Print statistics
        sum.reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: roi,
          scale: 30,
          maxPixels: 1e13
        }).evaluate(function(result) {
          print('🔹 Total water pixels (sum) for year ' + year + ':', result);
        });

        count.reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: roi,
          scale: 30,
          maxPixels: 1e13
        }).evaluate(function(result) {
          print('🔹 Total valid pixels (count) for year ' + year + ':', result);
        });

        fp.reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: roi,
          scale: 30,
          maxPixels: 1e13
        }).evaluate(function(result) {
          print('🔹 Water presence frequency (fp) for year ' + year + ':', result);
        });

        // Generate 2x2 patches for analysis
        var kernel = ee.Kernel.fixed(2, 2, [[1, 1], [1, 1]]);
        var arrayImage = binaryWater.neighborhoodToArray(kernel);
        
        var samples = arrayImage.sample({
          region: roi,
          scale: 30,
          numPixels: 100,
          seed: 42,
          geometries: true
        });

        print('🔍 2x2 patch samples:', samples.limit(5));

        // Export patch samples
        Export.table.toDrive({
          collection: samples,
          description: 'water_patches_' + year,
          fileFormat: 'CSV'
        });

      } else {
        print('⚠️ Invalid Otsu threshold for year ' + year);
      }
    });

    // Generate and display histogram analysis
    var histogramDict = mndwi.reduceRegion({
      reducer: ee.Reducer.histogram({ maxBuckets: 256 }),
      geometry: roi,
      scale: 30,
      bestEffort: true
    });

    var bandName = mndwi.bandNames().get(0);
    var histData = ee.Dictionary(histogramDict.get(bandName));
    var threshOtsu = getOtsuThreshold(mndwi, roi);

    threshOtsu.evaluate(function(thresh) {
      print('Extracted histogram for ' + year + ':');
      print('Otsu threshold for ' + year + ':', thresh);

      var histogramList = ee.List(histData.get('histogram'));
      var bucketMeansList = ee.List(histData.get('bucketMeans'));

      var histogram2D = histogramList.map(function(val) {
        return ee.List([val]);
      });

      var histogramChart = ui.Chart.array.values({
        array: ee.Array(histogram2D),
        axis: 0,
        xLabels: bucketMeansList
      })
      .setChartType('ColumnChart')
      .setOptions({
        title: null,
        hAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
        vAxis: { title: null, textPosition: 'none', gridlines: { color: 'transparent' } },
        legend: { position: 'none' },
        backgroundColor: { fill: 'transparent' },
        colors: ['#888888'],
        lineWidth: 1
      });

      print(histogramChart);
      print(ui.Label('⬆️ Otsu threshold for ' + year + ' ≈ ' + thresh.toFixed(3)));

      plotOtsuSeparationChart(histData, thresh, year);

      // Calculate class probabilities
      var totalCount = ee.Number(histogramList.reduce(ee.Reducer.sum()));

      var classA_sum = ee.Number(
        ee.List.sequence(0, bucketMeansList.length().subtract(1))
          .map(function(i) {
            i = ee.Number(i);
            var mean = ee.Number(bucketMeansList.get(i));
            var count = ee.Number(histogramList.get(i));
            return mean.lte(thresh) ? count : 0;
          })
          .reduce(ee.Reducer.sum())
      );

      var classB_sum = totalCount.subtract(classA_sum);
      var probA = classA_sum.divide(totalCount);
      var probB = classB_sum.divide(totalCount);

      plotOtsuProbabilities(probA, probB, year);
    });
  }
}

// Event listeners for UI controls
yearDropdown.onChange(updateMap);
checkboxRGB.onChange(updateMap);
checkboxNDWI.onChange(updateMap);
checkboxMNDWI.onChange(updateMap);
checkboxWater.onChange(updateMap);
checkboxLand.onChange(updateMap);
checkboxCoastline.onChange(updateMap);
checkboxOtsu.onChange(updateMap);

// ============================================================================
// SECTION 8: DATA EXPORT FUNCTIONS
// ============================================================================

// Utility functions for data export
function getDateRange(year) {
  return {
    startDate: ee.Date.fromYMD(year, 1, 1),
    endDate: ee.Date.fromYMD(year + 4, 12, 31)
  };
}

function exportImage(image, description) {
  Export.image.toDrive({
    image: image,
    description: description,
    folder: 'extract_LC',
    fileNamePrefix: description,
    scale: 30,
    region: roi,
    crs: 'EPSG:31983',
    maxPixels: 1e8
  });
}

// Tide data generation and export
function generateSpecificTideDates(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var highTideTimes = [
    [2, '13:00'], [2, '14:00'], [15, '13:00'],
    [15, '14:00'], [23, '13:00'], [23, '14:00'],
    [1, '13:00'], [1, '14:00'], [12, '13:00'],
    [12, '14:00'], [20, '13:00'], [20, '14:00']
  ];

  var dates = ee.List.sequence(0, 11).map(function(i) {
    var monthStartDate = startDate.advance(i, 'month');
    return ee.List(highTideTimes).map(function(dayTime) {
      var day = ee.Number(ee.List(dayTime).get(0)).subtract(1);
      var time = ee.List(dayTime).get(1);
      var date = monthStartDate.advance(day, 'day');
      var dateTime = date.format('YYYY-MM-dd').cat('T').cat(time);
      return ee.Feature(null, { 'date': dateTime });
    });
  }).flatten();

  return ee.FeatureCollection(dates);
}

// Metadata extraction functions
function exportMetadata(year) {
  var dateRange = getDateRange(year);
  var col = (year < 2013) ? l4.merge(l5).merge(l7) : l8.merge(l9);
  
  return col.filterBounds(roi)
            .filterDate(dateRange.startDate, dateRange.endDate);
}

var extractMetadata = function(image) {
  var metadata = ee.Image(image).toDictionary([
    'ALGORITHM_SOURCE_SURFACE_REFLECTANCE', 'ALGORITHM_SOURCE_SURFACE_TEMPERATURE', 'CLOUD_COVER', 'CLOUD_COVER_LAND', 
    'COLLECTION_CATEGORY', 'COLLECTION_NUMBER', 'CORRECTION_BIAS_BAND_1', 'CORRECTION_BIAS_BAND_2', 'CORRECTION_BIAS_BAND_3', 
    'CORRECTION_BIAS_BAND_4', 'CORRECTION_BIAS_BAND_5', 'CORRECTION_BIAS_BAND_6', 'CORRECTION_BIAS_BAND_7', 'CORRECTION_GAIN_BAND_1', 
    'CORRECTION_GAIN_BAND_2', 'CORRECTION_GAIN_BAND_3', 'CORRECTION_GAIN_BAND_4', 'CORRECTION_GAIN_BAND_5', 'CORRECTION_GAIN_BAND_6', 
    'CORRECTION_GAIN_BAND_7', 'DATA_SOURCE_AIR_TEMPERATURE', 'DATA_SOURCE_ELEVATION', 'DATA_SOURCE_OZONE', 'DATA_SOURCE_PRESSURE', 
    'DATA_SOURCE_REANALYSIS', 'DATA_SOURCE_WATER_VAPOR', 'DATA_TYPE_L0RP', 'DATE_ACQUIRED', 'DATE_PRODUCT_GENERATED', 'DATUM', 
    'EARTH_SUN_DISTANCE', 'ELLIPSOID', 'EPHEMERIS_TYPE', 'GEOMETRIC_RMSE_MODEL', 'GEOMETRIC_RMSE_MODEL_X', 'GEOMETRIC_RMSE_MODEL_Y', 
    'GEOMETRIC_RMSE_VERIFY', 'GEOMETRIC_RMSE_VERIFY_QUAD_LL', 'GEOMETRIC_RMSE_VERIFY_QUAD_LR', 'GEOMETRIC_RMSE_VERIFY_QUAD_UL', 
    'GEOMETRIC_RMSE_VERIFY_QUAD_UR', 'GRID_CELL_SIZE_REFLECTIVE', 'GRID_CELL_SIZE_THERMAL', 'GROUND_CONTROL_POINTS_MODEL', 
    'GROUND_CONTROL_POINTS_VERIFY', 'GROUND_CONTROL_POINTS_VERSION', 'IMAGE_QUALITY', 'L1_DATE_PRODUCT_GENERATED', 'L1_LANDSAT_PRODUCT_ID', 
    'L1_PROCESSING_LEVEL', 'L1_PROCESSING_SOFTWARE_VERSION', 'L1_REQUEST_ID', 'LANDSAT_PRODUCT_ID', 'LANDSAT_SCENE_ID', 'MAP_PROJECTION', 
    'ORIENTATION', 'PROCESSING_LEVEL', 'PROCESSING_SOFTWARE_VERSION', 'REFLECTANCE_ADD_BAND_1', 'REFLECTANCE_ADD_BAND_2', 
    'REFLECTANCE_ADD_BAND_3', 'REFLECTANCE_ADD_BAND_4', 'REFLECTANCE_ADD_BAND_5', 'REFLECTANCE_ADD_BAND_7', 'REFLECTANCE_MULT_BAND_1', 
    'REFLECTANCE_MULT_BAND_2', 'REFLECTANCE_MULT_BAND_3', 'REFLECTANCE_MULT_BAND_4', 'REFLECTANCE_MULT_BAND_5', 'REFLECTANCE_MULT_BAND_7', 
    'REFLECTIVE_LINES', 'REFLECTIVE_SAMPLES', 'REQUEST_ID', 'SATURATION_BAND_1', 'SATURATION_BAND_2', 'SATURATION_BAND_3', 
    'SATURATION_BAND_4', 'SATURATION_BAND_5', 'SATURATION_BAND_6', 'SATURATION_BAND_7', 'SCENE_CENTER_TIME', 'SENSOR_ANOMALIES', 
    'SENSOR_ID', 'SENSOR_MODE', 'SENSOR_MODE_SLC', 'SPACECRAFT_ID', 'STATION_ID', 'SUN_AZIMUTH', 'SUN_ELEVATION', 'TEMPERATURE_ADD_BAND_ST_B6', 
    'TEMPERATURE_MAXIMUM_BAND_ST_B6', 'TEMPERATURE_MINIMUM_BAND_ST_B6', 'TEMPERATURE_MULT_BAND_ST_B6', 'THERMAL_LINES', 'THERMAL_SAMPLES',
    'system:index', 'system:time_start', 'system:time_end', 'system:footprint', 'system:bands'
  ]);

  return ee.Feature(null, metadata);
};

// ============================================================================
// SECTION 9: EXPORT BUTTONS
// ============================================================================

// Landsat composite download
var downloadLandsatButton = ui.Button({
  label: 'Download Landsat Image',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);
    
    exportImage(image, 'Landsat_' + year);
  }
});

// NDWI download
var downloadNDWIButton = ui.Button({
  label: 'Download NDWI',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);

    var ndwi = image.expression('(B3 - B5) / (B3 + B5)', {
      B3: image.select('B3'),
      B5: image.select('B5')
    }).rename('NDWI');

    var visParams = {
      min: -1,
      max: 1,
      palette: ['brown', 'white', 'blue']
    };

    var ndwiRGB = ndwi.visualize(visParams);
    exportImage(ndwiRGB, 'NDWI_' + year);
  }
});

// MNDWI download
var downloadMNDWIButton = ui.Button({
  label: 'Download MNDWI',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);

    var mndwi = image.expression('(B3 - B5) / (B3 + B5)', {
      B3: image.select('B3'),
      B5: image.select('B5')
    }).rename('MNDWI');

    var visParams = {
      min: -1,
      max: 1,
      palette: ['blue', 'white', 'green']
    };

    var mndwiRGB = mndwi.visualize(visParams);
    exportImage(mndwiRGB, 'MNDWI_' + year);
  }
});

// Binary water (Otsu) download
var downloadOtsuButton = ui.Button({
  label: 'Download Binary Water (Otsu)',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);

    var mndwi = image.expression('(B3 - B5) / (B3 + B5)', {
      B3: image.select('B3'),
      B5: image.select('B5')
    }).rename('MNDWI');

    var otsuThreshold = getOtsuThreshold(mndwi, roi);

    otsuThreshold.evaluate(function(thresh) {
      if (thresh !== null) {
        var binaryWater = mndwi.gte(thresh).selfMask();
        exportImage(binaryWater, 'Binary_Water_Otsu_' + year);
      } else {
        print('Invalid or not found Otsu threshold.');
      }
    });
  }
});

// Water mask download
var downloadWaterButton = ui.Button({
  label: 'Download Water',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);

    var mndwi = image.expression('(B3 - B5) / (B3 + B5)', {
      B3: image.select('B3'),
      B5: image.select('B5')
    }).rename('MNDWI');
    
    var water = mndwi.gte(0.1).selfMask();
    exportImage(water, 'Water_' + year);
  }
});

// Land mask download
var downloadLandButton = ui.Button({
  label: 'Download Land',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);

    var mndwi = image.expression('(B3 - B5) / (B3 + B5)', {
      B3: image.select('B3'),
      B5: image.select('B5')
    }).rename('MNDWI');
    
    var land = mndwi.lt(0.1).selfMask();
    exportImage(land, 'Land_' + year);
  }
});

// Coastline download
var downloadCoastlineButton = ui.Button({
  label: 'Download Coastline',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var dateRange = getDateRange(year);
    var image = landsatComposite(roi, [dateRange.startDate, dateRange.endDate], year);

    var mndwi = image.expression('(B3 - B5) / (B3 + B5)', {
      B3: image.select('B3'),
      B5: image.select('B5')
    }).rename('MNDWI');
    
    var coastline = ee.Algorithms.CannyEdgeDetector({
      image: mndwi, threshold: 0.7, sigma: 1
    }).selfMask();
    
    exportImage(coastline, 'Coastline_' + year);
  }
});

// Metadata download
var downloadMetadataButton = ui.Button({
  label: 'Download Landsat Metadata',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var filteredCollection = exportMetadata(year);
    var metadata = filteredCollection.map(extractMetadata);

    Export.table.toDrive({
      collection: metadata,
      description: 'Landsat_Metadata_' + year,
      folder: 'extract_LC',
      fileFormat: 'CSV'
    });
  }
});

// Tide data download
var downloadTideDataButton = ui.Button({
  label: 'Download Tide Data',
  onClick: function() {
    var selectedYear = yearDropdown.getValue();
    if (!selectedYear) return;

    var year = parseInt(selectedYear, 10);
    var tideDates = generateSpecificTideDates(year);

    Export.table.toDrive({
      collection: tideDates,
      description: 'Tide_Dates_' + year,
      folder: 'extract_LC',
      fileFormat: 'CSV'
    });
  }
});

// ============================================================================
// SECTION 10: FINAL UI ASSEMBLY
// ============================================================================

// Add controls to panel
panel.add(ui.Label('Choose the year:'));
panel.add(yearDropdown);

// Add checkboxes
panel.add(ui.Label('Layer Visibility:'));
panel.add(checkboxRGB);
panel.add(checkboxNDWI);
panel.add(checkboxMNDWI);
panel.add(checkboxOtsu);
panel.add(checkboxWater);
panel.add(checkboxLand);
panel.add(checkboxCoastline);

// Add download buttons
panel.add(ui.Label('Select products to download:'));
panel.add(downloadLandsatButton);
panel.add(downloadNDWIButton);
panel.add(downloadMNDWIButton);
panel.add(downloadOtsuButton);
panel.add(downloadWaterButton);
panel.add(downloadLandButton);
panel.add(downloadCoastlineButton);
panel.add(downloadMetadataButton);
panel.add(downloadTideDataButton);

// Add footer
var footer = ui.Label({
  value: 'Click RUN to execute the code. Use the buttons above to select and download products.',
  style: {fontSize: '14px', color: 'gray', textAlign: 'center', margin: '10px'}
});
panel.add(footer);

// Add panel to UI
ui.root.insert(0, panel);

// Final status message
print('🚀 Landsat Shoreline Analysis Interface loaded successfully!');
print('📊 All products are ready for visualization and export.');

// Citation reminder
// alert('When citing, use: "DIAS, K. S. & MENEZES, R. A. A., 2025."');
