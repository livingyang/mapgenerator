(function () {

function Point(x, y) {
    this.x = x;
    this.y = y;  
};

function Line(pointA, pointB) {
    this.points = new Array();
    
    this.points.push(pointA);
    this.points.push(pointB);
};

function Hexagon(lineA, lineB, lineC, lineD, lineE, lineF) {
    this.used = false;
    this.lines = new Array();
    this.outline = new Array();
    this.neighbors = new Array();
    
    this.lines.push(lineA);
    this.lines.push(lineB);
    this.lines.push(lineC);
    this.lines.push(lineD);
    this.lines.push(lineE);
    this.lines.push(lineF);
};

/* Helpers and native extentions */
function distance(a, b) {
    return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
};

Array.prototype.extend = function (array) {
    for (var i = 0; i < array.length; i++) {
        this.push(array[i]);
    }
    
    return this;
};

if (!Array.prototype.indexOf) {
    Array.prototype.indexOf = function(item, from) {
        var len = this.length;
        
        for (var i = (from < 0) ? Math.max(0, len + from) : from || 0; i < len; i++) {
            if (this[i] === item)
                return i;
        }
        
        return -1;   
    };
}

Array.prototype.contains = function(item, from) {
    return this.indexOf(item, from) != -1;
};

Array.prototype.include = function(item) {
    if (!this.contains(item))
        this.push(item);
    
    return this;  
};

Array.prototype.erase = function(item) {
    for (var i = this.length; i--; i) {
        if (this[i] === item)
            this.splice(i, 1);
    }
    
    return this;  
};

Array.prototype.getLast = function() {
    return (this.length) ? this[this.length - 1] : null;
};

function Country() {
    this.hexagons = new Array();
    this.neighbors = new Array();
    this.outline = new Array();
    this.inlines = new Array();
};

Country.prototype.getRandomNeighborHexagon = function(useCompactShapes) {
    if (useCompactShapes) {
        while (true) {
            var hexagon = this.hexagons[rand(0, this.hexagons.length - 1)];
            var neighborHexagon = hexagon.neighbors[rand(0, hexagon.neighbors.length - 1)];
                
            if (!neighborHexagon.used)
                return neighborHexagon;
        }
    }
    else {
        var neighborHexagons = new Array();
        
        for (var i = 0; i < this.hexagons.length; i++) {
            for (var j = 0; j < this.hexagons[i].neighbors.length; j++) {
                if (!this.hexagons[i].neighbors[j].used)
                    neighborHexagons.include(this.hexagons[i].neighbors[j]);
            }
        }
        
        return neighborHexagons[rand(0, neighborHexagons.length -1)];
    }
};

Country.prototype.getPointField = function(points, lineLength) {
    var connectedPoints = new Array(),
        isConnected = true;
    
    connectedPoints.push(points[0]);
    points.erase(points[0]);
    
    while (isConnected) {
        isConnected = false;
        
        for (var i = 0; i < connectedPoints.length; i++) {
            for (var j = 0; j < points.length; j++) {
                if (distance(connectedPoints[i], points[j]) < lineLength) {
                    connectedPoints.push(points[j]);
                    points.erase(points[j]);
                    
                    isConnected = true;
                    break;
                }
            }
            
            if (isConnected)
                break;
        }
    }
    
    return connectedPoints;
};

Country.prototype.getHexagonField = function(hexagons) {
    var connectedHexagons = new Array(),
        isConnected = true;
    
    connectedHexagons.push(hexagons[0]);
    hexagons.erase(hexagons[0]);
    
    while (isConnected) {
        isConnected = false;
        
        for (var i = 0, ii = connectedHexagons.length; i < ii; i++) {
            for (var j = 0, jj = hexagons.length; j < jj; j++) {
                if (hexagons[j].neighbors.contains(connectedHexagons[i])) {
                    connectedHexagons.push(hexagons[j]);
                    hexagons.erase(hexagons[j]);
                    
                    isConnected = true;
                    break;
                }
            }
            
            if (isConnected)
                break;
        }
    }
    
    return connectedHexagons;
};

Country.prototype.getCenter = function() {
    // triplePoints are points in the inside of a country / a triplePoint is part of 3 hexagons
    var triplePoints = new Array(),
        points = new Array();
    
    for (var i = 0, ii = this.inlines.length; i < ii; i++) {
        for (var j = 0; j < 2; j++) {
            var point = this.inlines[i].points[j];
            
            if (points.contains(point)) {
                if (!triplePoints.contains(point))
                    triplePoints.push(point);
            }
            else
                points.push(point); 
        }
    }
    
    var sumX = 0,
        sumY = 0,
        line = this.hexagons[0].lines[0];
        lineLength = distance(line.points[0], line.points[1]) * 1.1;
    
    // no triplePoints: set center in middle of a hexagon
    if (triplePoints.length < 1) {
        for (var i = 0, ii = this.hexagons.length; i < ii; i++) {
            var inCountryNeighbors = 0;
            
            for (var j = 0; j < 6; j++) {
                if (this.hexagons.contains(this.hexagons[i].neighbors[j])) {
                    inCountryNeighbors++;
                    
                    // use a hexagon as center with 3 neighbors
                    if (inCountryNeighbors == 3) {
                        for (var k = 0; k < 6; k++) {
                            sumX += this.hexagons[i].lines[k].points[0].x + this.hexagons[i].lines[k].points[1].x;
                            sumY += this.hexagons[i].lines[k].points[0].y + this.hexagons[i].lines[k].points[1].y;
                        }
                        
                        this.center = new Point(sumX / 12, sumY / 12);
                        return;
                    }
                }
            }
        }
        
        // if there is no hexagon with 3 neighbors, use first
        for (var i = 0; i < 6; i++) {
            sumX += this.hexagons[0].lines[i].points[0].x + this.hexagons[0].lines[i].points[1].x;
            sumY += this.hexagons[0].lines[i].points[0].y + this.hexagons[0].lines[i].points[1].y;
        }
        
        this.center = new Point(sumX / 12, sumY / 12);
        return;
    }
    
    // a pointField consists of connected triplePoints
    var pointFields = new Array();
    
    while (triplePoints.length > 0) {
        pointFields.push(this.getPointField(triplePoints, lineLength));
    }
    
    // get the biggest pointField
    var pointField = pointFields[0];
    for (var i = 1; i < pointFields.length; i++) {
        if (pointFields[i].length > pointField.length)
            pointField = pointFields[i];
    }
    
    // inLineHexagons are hexagons completely on the inside of a country
    var inlineHexagons = new Array();
    
    for (var i = 0, ii = this.hexagons.length; i < ii; i++) {
        var containsHex = true;
        
        for (var j = 0; j < 6; j++) {
            if (!pointField.contains(this.hexagons[i].lines[j].points[0]) || 
                !pointField.contains(this.hexagons[i].lines[j].points[1])) {
                
                containsHex = false;
                break;
            }
        }
        
        if (containsHex) 
            inlineHexagons.push(this.hexagons[i]);
    }
    
    // no inlineHexagons: set center to average triplePoint in pointField
    if (inlineHexagons.length < 1) {
        for (var i = 0, ii = pointField.length; i < ii; i++) {
            sumX += pointField[i].x;
            sumY += pointField[i].y;
        }
        
        var centerPoint = new Point(sumX / pointField.length, sumY / pointField.length),
            minDistance = Infinity,
            j;
        
        if (pointField.length < 7) {
            this.center = centerPoint;
            return;
        }
        
        for (var i = 0, ii = pointField.length; i < ii; i++) {
            var pointDistance = distance(centerPoint, pointField[i]);
            
            if (pointDistance < minDistance) {
                j = i;
                minDistance = pointDistance;
                
                if (minDistance < lineLength / 2)
                    break;
            }
        }
        
        this.center = pointField[j];
        return;
    }
    
    // a hexagonField is a field consisting of inlineHexagons
    var hexagonFields = new Array();
    
    while (inlineHexagons.length > 0) {
        hexagonFields.push(this.getHexagonField(inlineHexagons));
    }
    
    // get the biggest hexagonField
    var hexagonField = hexagonFields[0];
    
    for (var i = 1; i < hexagonFields.length; i++) {
        if (hexagonFields[i].length > hexagonField.length)
            hexagonField = hexagonFields[i];
    }
    
    var hexagonCenters = new Array();
    
    for (var i = 0, ii = hexagonField.length; i < ii; i++) {
        var x = 0,
            y = 0;
            
        for (var j = 0; j < 6; j++) {    
            x += hexagonField[i].lines[j].points[0].x + hexagonField[i].lines[j].points[1].x,
            y += hexagonField[i].lines[j].points[0].y + hexagonField[i].lines[j].points[1].y;
        }
        
        sumX += x /= 12;
        sumY += y /= 12;
        hexagonCenters.push(new Point(x, y));
    }
    
    var centerPoint = new Point(sumX / hexagonField.length, sumY / hexagonField.length),
        minDistance = Infinity,
        j;
    
    if (hexagonField.length < 7) {
        this.center = centerPoint;
        return;
    }
    
    for (var i = 0, ii = hexagonField.length; i < ii; i++) {
        var hexDistance = distance(centerPoint, hexagonCenters[i]);
        
        if (hexDistance < minDistance) {
            j = i;
            minDistance = hexDistance;
            
            if (minDistance < lineLength * 2 / 3)
                break;
        }
    }
    
    this.center = new Point(hexagonCenters[j].x, hexagonCenters[j].y);
    return;
};

Country.prototype.generateOutline = function() {
    // lineArray containing only outlines
    var outLines = new Array();
    
    for (var i = 0, ii = this.hexagons.length; i < ii; i++) {
        for (var j = 0; j < 6; j++) {
            var line = this.hexagons[i].lines[j];
            
            if (outLines.contains(line)) {
                outLines = outLines.erase(line);
                this.inlines.push(line);
            }
            else
                outLines.push(line);
        }
    }
    
    // getting line on the outside
    var line = outLines.getLast();
    for (var i = 0; i < outLines.length; i++) {
        if (outLines[i].points[0].x < line.points[0].x)
            line = outLines[i];
    }
    
    // creating the outline and bounding box
    this.outline.push(line.points[0]);
    this.outline.push(line.points[1]);
    outLines = outLines.erase(line);
    
    var startPoint = line.points[0],
        point = line.points[1];
    
    this.boundingBox = {};
    this.boundingBox.min = new Point(Math.min(startPoint.x, point.x), Math.min(startPoint.y, point.y));
    this.boundingBox.max = new Point(Math.max(startPoint.x, point.x), Math.max(startPoint.y, point.y));
    
    while (startPoint != point) {
        for (var i = 0; i < outLines.length; i++) {
            var isNewPoint = false;
            
            if (outLines[i].points[0] == point) {   
                point = outLines[i].points[1];
                this.outline.push(outLines[i].points[1]);
                outLines = outLines.erase(outLines[i]);
                isNewPoint = true;
            } 
            else if (outLines[i].points[1] == point) {
                point = outLines[i].points[0];
                this.outline.push(outLines[i].points[0]);
                outLines = outLines.erase(outLines[i]);
                isNewPoint = true;
            }
            
            if (isNewPoint) {
                this.boundingBox.min.x = Math.min(this.boundingBox.min.x, point.x);
                this.boundingBox.min.y = Math.min(this.boundingBox.min.y, point.y);
                
                this.boundingBox.max.x = Math.max(this.boundingBox.max.x, point.x);
                this.boundingBox.max.y = Math.max(this.boundingBox.max.y, point.y);
            }
        }
    }
    
    if (outLines.length > 0)
        this.holeLines = outLines;
};

/* Helpers and native extentions */
function rand(minimum, maximum) {
    return Math.floor(Math.random() * (maximum - minimum + 1)) + minimum;
}

if (!Array.prototype.indexOf) {
    Array.prototype.indexOf = function(item, from) {
        var len = this.length;
        
        for (var i = (from < 0) ? Math.max(0, len + from) : from || 0; i < len; i++) {
            if (this[i] === item)
                return i;
        }
        
        return -1;   
    };
}

Array.prototype.contains = function(item, from) {
    return this.indexOf(item, from) != -1;
};

Array.prototype.include = function(item) {
    if (!this.contains(item))
        this.push(item);
    
    return this;  
};

Array.prototype.combine = function(array) {
    for (var i = 0, length = array.length; i < length; i++) {
        this.include(array[i]);
    }
    
    return this;
};

if (!Array.prototype.forEach) {
    Array.prototype.forEach = function(fun /*, thisp*/) {
        var len = this.length >>> 0;
        
        if (typeof fun != "function")
            throw new TypeError();
    
        var thisp = arguments[1];
        
        for (var i = 0; i < len; i++) {
            if (i in this)
                fun.call(thisp, this[i], i, this);
        }
    };
}

function Map(width, height, hexagonSize) {
    this.points = new Array();
    this.lines = new Array();
    this.hexagons = new Array();
    this.usedHexagons = new Array();
    this.countries = new Array();
    this.width = width;
    this.height = height;
    this.hexagonSize = hexagonSize;
};

Map.prototype.clear = function() {
    if (this.countries.length > 0) {
        this.usedHexagons = new Array();
        this.countries = new Array();
    
        for (var i = 0, ii = this.hexagons.length; i < ii; i++) {
            this.hexagons[i].used = false;
        }
    }
};

Map.prototype.generateHexagonArray = function(useDistortion) {
    var hexagonWidth = Math.sqrt(3) * this.hexagonSize / 2,
        numberOfHexagonsInARow = parseInt((this.width / hexagonWidth) - 0.5), 
        numberOfHexagonsInAColumn = parseInt(((4 * this.height) / (3 * this.hexagonSize)) - (1 / 3)),
        distortionAmount = 1;
    
    // pointArray
    for (var row = 0; row < numberOfHexagonsInAColumn + 1; row++) {
        for (var column = 0; column < numberOfHexagonsInARow + 1; column++) {
            var x, y, phi, r;
            
            x = column * hexagonWidth;
            y = row * this.hexagonSize * 0.75;
            
            if ((row % 2) == 0)
                y += 0.25 * this.hexagonSize;
             
            if (useDistortion) {   
                phi = Math.random() * Math.PI * 2;
                r = Math.random() * this.hexagonSize/4 * distortionAmount;
                x += r * Math.cos(phi);
                y += r * Math.sin(phi);
            }
            
            this.points.push(new Point(x, y));
                
            x = (column + 0.5) * hexagonWidth;
            y = row * this.hexagonSize * 0.75;
            
            if ((row % 2) == 1)
                y += 0.25 * this.hexagonSize;
            
            if (useDistortion) {   
                phi = Math.random() * 2 * Math.PI;
                r = Math.random() * this.hexagonSize/4 * distortionAmount;
                x += r * Math.cos(phi);
                y += r * Math.sin(phi);
            }
                
            this.points.push(new Point(x, y));
        }
    }
    
    // lineArray
    for (var row = 0; row < numberOfHexagonsInAColumn + 1; row++) {
        var number = numberOfHexagonsInARow * 2 + 2;
        
        for (var column = 0; column < numberOfHexagonsInARow * 2 + 1; column++) {
            var pointA = this.points[(row * number) + column],
                pointB = this.points[(row * number) + column + 1];
            this.lines.push(new Line(pointA, pointB));
        }
        
        if (row < numberOfHexagonsInAColumn) {
            var oddrow = (row % 2) ? 1 : 0;
            
            for (var column = 0; column < numberOfHexagonsInARow + 1; column++) {
                var pointA = this.points[(row * number) + column * 2 + oddrow],
                    pointB = this.points[((row + 1) * number) + column * 2 + oddrow];
                this.lines.push(new Line(pointA, pointB));
            }
        }
    }
    
    // hexagonArray
    var linesPerRow = numberOfHexagonsInARow * 3 + 2;
    
    for (var row = 0; row < numberOfHexagonsInAColumn; row++) {
        var oddrow = (row % 2) ? 1 : 0;
            
        for (var column = 0; column < numberOfHexagonsInARow; column++) {
            var number = linesPerRow * row + column * 2 + oddrow;
            var lineA = this.lines[number],
                lineB = this.lines[number+1];
            
            number = linesPerRow * row + (numberOfHexagonsInARow * 2 + 1) + column;
            
            var lineC = this.lines[number],
                lineD = this.lines[number+1];
            
            number = linesPerRow * (row + 1) + column * 2 + oddrow;
            
            var lineE = this.lines[number],
                lineF = this.lines[number+1];
            
            this.hexagons.push(new Hexagon(lineA, lineB, lineD, lineF, lineE, lineC));
        }
    }
    
    // hexagonNeighbors
    var leftBorder = true,
        topBorder = true,
        rightBorder = false,
        bottomBorder = false,
        index = 0;
    
    for (var i = 0; i < numberOfHexagonsInAColumn; i++) {
        leftBorder = true;
        
        if (i == (numberOfHexagonsInAColumn - 1))
            bottomBorder = true;
        
        for (var j = 0; j < numberOfHexagonsInARow; j++) {
            if (j == (numberOfHexagonsInARow - 1))
                rightBorder = true;
            
            if (!leftBorder)
                this.hexagons[index].neighbors.push(this.hexagons[index - 1]);
            
            if (!rightBorder)
                this.hexagons[index].neighbors.push(this.hexagons[index + 1]);
                
            if (!topBorder) {
                this.hexagons[index].neighbors.push(this.hexagons[index - numberOfHexagonsInARow]);
                
                if ((i % 2) == 1) {
                    if (!rightBorder)
                        this.hexagons[index].neighbors.push(this.hexagons[index + 1 - numberOfHexagonsInARow]);
                } else {
                    if (!leftBorder)
                        this.hexagons[index].neighbors.push(this.hexagons[index - 1 - numberOfHexagonsInARow]);
                }
            }
            
            if (!bottomBorder) {
                this.hexagons[index].neighbors.push(this.hexagons[index + numberOfHexagonsInARow]);
                
                if ((i % 2) == 1) {
                    if (!rightBorder)
                        this.hexagons[index].neighbors.push(this.hexagons[index + 1 + numberOfHexagonsInARow]);
                } else {
                    if (!leftBorder)
                        this.hexagons[index].neighbors.push(this.hexagons[index - 1 + numberOfHexagonsInARow]);
                }
            }
                
            if (leftBorder)
                leftBorder = false;
            else if (rightBorder)
                rightBorder = false;
                
            index++;
        }
        
        if (topBorder)
            topBorder = false;
    }
};

Map.prototype.holeChecker = function(hexagon, maximumHoleSize) {
    var freeHexagons = new Array();
    
    freeHexagons.push(hexagon);
    
    for (var i = 0; i < freeHexagons.length; i++) {
        if (freeHexagons.length >= maximumHoleSize)
            return false;
        
        for (var j = 0; j < freeHexagons[i].neighbors.length; j++) {
            if (!freeHexagons[i].neighbors[j].used) {
                freeHexagons.include(freeHexagons[i].neighbors[j]);
            }
        }
    }
    
    if (freeHexagons.length >= maximumHoleSize)
        return false;
    else {
        this.usedHexagons.combine(freeHexagons);
        freeHexagons.forEach(function(hexagon) {
            hexagon.used = true;
        });
        
        return true;
    }    
};

Map.prototype.generateCountry = function(ID, neighborCountry, size, useCompactShapes) {
    var country = new Country(),
        startHexagon;
    
    country.ID = ID;
    
    if (neighborCountry != null) {
        do {
            startHexagon = neighborCountry.getRandomNeighborHexagon(true);
            
            if (!startHexagon)
                throw 'Epic Fail';
            
        } while(this.holeChecker(startHexagon, size))
    }
    else
        startHexagon = this.hexagons[rand(0, this.hexagons.length - 1)];
      
    startHexagon.used = true;    
    country.hexagons.push(startHexagon);
    this.usedHexagons.push(startHexagon);
    
    for (var i = 1; i < size; i++) {
        var newHexagon = country.getRandomNeighborHexagon(useCompactShapes);
        
        newHexagon.used = true;
        country.hexagons.push(newHexagon);
        this.usedHexagons.push(newHexagon);
    }
    
    return country;
};

Map.prototype.normalGenerator = function(numberOfCountries, countrySizeVariance, useCompactShapes) {
    var mapCoverage = 0.6;
    
    var averageCountrySize = parseInt((this.hexagons.length - this.usedHexagons.length) * mapCoverage / numberOfCountries);
    
    if (countrySizeVariance < 0 || countrySizeVariance > 0.9)
        countrySizeVariance = 0;
    
    for (var i = 0; i < numberOfCountries; i++) {
        var countrySize = (averageCountrySize + rand(0, parseInt(averageCountrySize * countrySizeVariance)) * (rand(0, 1) ? 1 : -1));
        
        if (this.countries.length > 0) {
            var globalCountry = new Country();
            
            globalCountry.hexagons = this.usedHexagons;
            this.countries.push(this.generateCountry(i, globalCountry, countrySize, useCompactShapes));
        }
        else 
            this.countries.push(this.generateCountry(i, null, countrySize, useCompactShapes));   
    }
};

Map.prototype.getCountryNeighbors = function() {
    for (var i = 0, ii = this.countries.length; i < ii; i++) {
        var countryOutline = this.countries[i].outline;
        
        for (var j = i + 1; j < ii; j++) {
            var minXIndex = this.countries[i].boundingBox.min.x < this.countries[j].boundingBox.min.x ? i : j,
                minYIndex = this.countries[i].boundingBox.min.y < this.countries[j].boundingBox.min.y ? i : j;
            
            if (this.countries[minXIndex].boundingBox.max.x >= this.countries[minXIndex == j ? i : j].boundingBox.min.x &&
                this.countries[minYIndex].boundingBox.max.y > this.countries[minYIndex == j ? i : j].boundingBox.min.y) {
                
                for (var k = 0, kk = this.countries[j].outline.length; k < kk; k += 2) {
                    if (countryOutline.contains(this.countries[j].outline[k])) {
                        this.countries[i].neighbors.push(this.countries[j]);
                        this.countries[j].neighbors.push(this.countries[i]);
                        break;
                    }
                }
            }
        }
    } 
};

Map.prototype.deleteCountryHoles = function() {
    for (var i = 0, ii = this.countries.length; i < ii; i++) {
        if (this.countries[i].holeLines != undefined) {
            var country = this.countries[i],
                holeHexagons = new Array();
            
            while (country.holeLines.length > 0) {
                for (var j = 0, jj = this.hexagons.length; j < jj; j++) {
                    if (this.hexagons[j].lines.contains(country.holeLines[0]) && 
                        !country.hexagons.contains(this.hexagons[j])) {
                        holeHexagons.push(this.hexagons[j]);
                        break;
                    }
                }
                
                while (holeHexagons.length > 0) {
                    var hexagon = holeHexagons.pop();
                    country.hexagons.push(hexagon);
                    
                    for (var j = 0; j < 6; j++) {
                        country.inlines.include(hexagon.lines[j]);
                        country.holeLines.erase(hexagon.lines[j]);
                    }
                    
                    for (var j = 0, jj = hexagon.neighbors.length; j < jj; j++) {
                        if (!country.hexagons.contains(hexagon.neighbors[j]))
                            holeHexagons.push(hexagon.neighbors[j]);
                    }
                }
            }
        }
    }
};

Map.prototype.calculateOutlines = function() {
    for (var i = 0; i < this.countries.length; i++) {
        this.countries[i].generateOutline();
    }        
};
    
Map.prototype.calculateCenters = function() {
    for (var i = 0; i < this.countries.length; i++) {
        this.countries[i].getCenter();
    }
};

function MapGenerator() {
};

MapGenerator.prototype.createHexagonPattern = function(mapWidth, mapHeight, hexagonSize, useDistortion) {
    this.map = new Map(mapWidth, mapHeight, hexagonSize);
    this.map.generateHexagonArray(useDistortion);
};

MapGenerator.prototype.generate = function(numberOfCountries, countrySizeVariance, useCompactShapes) {
    if (this.map == undefined)
        throw "call MapGenerator.createHexagonPattern() before generating";
    
    this.map.clear();
    this.map.normalGenerator(numberOfCountries, countrySizeVariance, useCompactShapes);
    this.map.calculateOutlines();
    this.map.deleteCountryHoles();
    this.map.calculateCenters();
    this.map.getCountryNeighbors();
};
    
MapGenerator.prototype.getCountries = function() {
    return this.map.countries;
};

MapGenerator.prototype.getRawMap = function() {
    return this.map;  
};

MapGenerator.prototype.getMap = function() {
    var map = {};
    
    map.width = this.map.width;
    map.height = this.map.height;
    map.regions = new Array();
    
    for (var i = 0; i < this.map.countries.length; i++) {
        var region = {};
        
        region.center = this.map.countries[i].center;
        region.ID = this.map.countries[i].ID;
        
        var pathString = "M " + this.map.countries[i].outline[0].x + " " + this.map.countries[i].outline[0].y;
        
        for (var j = 1; j < this.map.countries[i].outline.length; j++) {
            pathString += "L " + this.map.countries[i].outline[j].x + " " + this.map.countries[i].outline[j].y;
        }
        
        pathString += " Z";
        region.pathString = pathString;
        
        region.neighborIDs = new Array();
        for (var j = 0; j < this.map.countries[i].neighbors.length; j++) {
            region.neighborIDs.push(this.map.countries[i].neighbors[j].ID)
        }
        
        map.regions.push(region);
    }
    
    map.adjacencyMatrix = new Array();
    
    for (var i = 0; i < this.map.countries.length; i++) {
        map.adjacencyMatrix[i] = new Array();
    
        for (var j = 0; j < this.map.countries.length; j++) {
            map.adjacencyMatrix[i][j] = 0;
        }
    }
    
    for (var i = 0; i < this.map.countries.length; i++) {
        for (var j = 0; j < this.map.countries[i].neighbors.length; j++) {
            var differenceX = this.map.countries[i].center.x - this.map.countries[i].neighbors[j].center.x;
            var differenceY = this.map.countries[i].center.y - this.map.countries[i].neighbors[j].center.y;
            var distance = Math.sqrt(Math.pow(differenceX, 2) + Math.pow(differenceY, 2));
        
            map.adjacencyMatrix[this.map.countries[i].ID][this.map.countries[i].neighbors[j].ID] = distance;
        }
    }
    
    return map;
};

this.MapGenerator = MapGenerator;

}).call(this);
