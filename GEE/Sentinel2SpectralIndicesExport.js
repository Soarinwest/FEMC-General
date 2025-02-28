/**
 * Script to process Sentinel-2 data for selected states,
 * calculate various spectral indices, and export the results.
 *
 * Requirements: Google Earth Engine account.
 * Author: Soren Donisvitch
 * Date: 01/15/2025
 * Foreword: The use or application of these code without permission of the author is prohibited.The author is not ->
 *  liable for the use, modification, or any other application of this or other provided scripts.
 */

// -------------------------------
// 1. Define Regions & Load Datasets
// -------------------------------

// List of states to process
var stateNames = ['Connecticut', 'Maine', 'Massachusetts', 'New Hampshire', 'Rhode Island', 'Vermont', 'New York'];

// Load the TIGER states dataset and filter for the selected states
var states = ee.FeatureCollection('TIGER/2018/States');
var northeastStates = states.filter(ee.Filter.inList('NAME', stateNames));

// Center the map on the selected state(s) and add a layer for visualization
Map.centerObject(northeastStates.geometry(), 6);
Map.addLayer(northeastStates, {color: 'blue'}, 'Selected States');

// -------------------------------
// 2. Define Helper Functions
// -------------------------------

/**
 * Applies a cloud mask to Sentinel-2 images using the QA60 band.
 * Bits 10 and 11 correspond to clouds and cirrus.
 */
function cloudMaskS2H(image) {
  var qa = image.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0)
              .and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask);
}

/**
 * Calculates several spectral indices and adds them as bands.
 */
function calculateIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = image.expression(
    '2.5 * ((B8 - B4) / (B8 + 6 * B4 - 7.5 * B2 + 1))',
    {
      'B8': image.select('B8'),
      'B4': image.select('B4'),
      'B2': image.select('B2')
    }
  ).rename('EVI');
  var savi = image.expression(
    '((B8 - B4) / (B8 + B4 + 0.5)) * 1.5',
    {
      'B8': image.select('B8'),
      'B4': image.select('B4')
    }
  ).rename('SAVI');
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var bsi = image.expression(
    '((B11 + B4) - (B8 + B2)) / ((B11 + B4) + (B8 + B2))',
    {
      'B11': image.select('B11'),
      'B4': image.select('B4'),
      'B8': image.select('B8'),
      'B2': image.select('B2')
    }
  ).rename('BSI');
  var ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI');
  var mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  var nbr = image.normalizedDifference(['B8', 'B12']).rename('NBR');

  return image.addBands([ndvi, evi, savi, ndwi, bsi, ndbi, mndwi, nbr]);
}

/**
 * Computes a rolling mean for an ImageCollection over a specified window size.
 */
function computeRollingMean(collection, windowSize) {
  var halfWindow = ee.Number(windowSize).divide(2);
  return collection.map(function(image) {
    var date = image.date();
    var start = date.advance(halfWindow.multiply(-1), 'day');
    var end = date.advance(halfWindow, 'day');
    var window = collection.filterDate(start, end);
    return window.mean().copyProperties(image, ['system:time_start']);
  });
}

// -------------------------------
// 3. Define Processing Parameters
// -------------------------------

// List of spectral indices for export
var spectralIndices = ['NDVI', 'EVI', 'SAVI', 'NDWI', 'BSI', 'NDBI', 'MNDWI', 'NBR'];

// Define the years and seasonal date range for analysis
var years = ee.List.sequence(2016, 2024);
var startMonth = 5; // May
var endMonth = 10;  // October

// -------------------------------
// 4. Main Processing Function
// -------------------------------

/**
 * Processes a given region by loading Sentinel-2 data, calculating indices,
 * applying temporal smoothing, and exporting each index as an image.
 *
 * @param {Geometry} regionGeometry - The region to process.
 * @param {string} regionName - Name to label the region.
 */
function processRegion(regionGeometry, regionName) {
  // Visualize the region on the map
  Map.addLayer(regionGeometry, {}, regionName);

  years.getInfo().forEach(function(y) {
    var year = ee.Number(y);
    var startDate = ee.Date.fromYMD(year, startMonth, 1);
    var endDate = ee.Date.fromYMD(year, endMonth, 31);

    // Load and filter Sentinel-2 data over the region and time period
    var s2Collection = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
      .filterBounds(regionGeometry)
      .filterDate(startDate, endDate)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
      .map(cloudMaskS2H)
      .map(calculateIndices);

    // Apply temporal smoothing using a 20-day rolling mean
    var smoothedCollection = computeRollingMean(s2Collection, 20);

    // Create a median composite of the spectral indices
    var indicesMedian = smoothedCollection.select(spectralIndices).median()
      .set('year', year)
      .toFloat()
      .clip(regionGeometry);

    // Export each spectral index separately to Google Drive
    spectralIndices.forEach(function(index) {
      Export.image.toDrive({
        image: indicesMedian.select(index).rename(index),
        description: 'Indices_' + index + '_' + regionName.replace(' ', '_') + '_' + year.getInfo(),
        folder: 'Northeast_Spectral_Indices',
        scale: 10,
        maxPixels: 1e13,
        region: regionGeometry,
        crs: 'EPSG:32618'
      });
    });
  });
}

// -------------------------------
// 5. Process Each State and Split Regions When Needed
// -------------------------------

stateNames.forEach(function(stateName) {
  var state = northeastStates.filter(ee.Filter.eq('NAME', stateName));
  var roi = state.geometry();

  // For states that need splitting (New York and Maine)
  if (stateName === 'New York' || stateName === 'Maine') {
    // Compute bounding box coordinates for splitting
    var bounds = roi.bounds();
    var coords = ee.List(bounds.coordinates().get(0));
    var xmin = ee.Number(ee.List(coords.get(0)).get(0));
    var ymin = ee.Number(ee.List(coords.get(0)).get(1));
    var xmax = ee.Number(ee.List(coords.get(2)).get(0));
    var ymax = ee.Number(ee.List(coords.get(2)).get(1));

    if (stateName === 'New York') {
      // Split New York into East and West regions along the central longitude
      var midLon = xmin.add(xmax).divide(2);
      var westRegion = roi.intersection(ee.Geometry.Rectangle([xmin, ymin, midLon, ymax]), ee.ErrorMargi (1));
      var eastRegion = roi.intersection(ee.Geometry.Rectangle([midLon, ymin, xmax, ymax]), ee.ErrorMargin(1));

      processRegion(westRegion, 'New_York_West');
      processRegion(eastRegion, 'New_York_East');
    } else if (stateName === 'Maine') {
      // Split Maine into North and South regions along the central latitude
      var midLat = ymin.add(ymax).divide(2);
      var southRegion = roi.intersection(ee.Geometry.Rectangle([xmin, ymin, xmax, midLat]), ee.ErrorMargin(1));
      var northRegion = roi.intersection(ee.Geometry.Rectangle([xmin, midLat, xmax, ymax]), ee.ErrorMargin(1));

      processRegion(southRegion, 'Maine_South');
      processRegion(northRegion, 'Maine_North');
    }
  } else {
    // Process the state as a single region if no split is needed
    processRegion(roi, stateName);
  }
});