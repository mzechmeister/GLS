var lines;

function getCSV(url, func) {
   var rawFile = new XMLHttpRequest();
 
   rawFile.open("GET", url, false);
   rawFile.overrideMimeType("text/plain");  // otherwise in firefox XML Parsing Error: syntax error
   rawFile.onreadystatechange = function () {
      if (rawFile.readyState === 4)
         if (rawFile.status === 200 || rawFile.status == 0)
            lines = rawFile.responseText;
            func();
        };
        rawFile.send();
}

function handleURL(url) {
   // Check for the various File API support.
   console.log("requesting", url)
   getCSV(url, function(){
      if (lines.length)  {
         lines = lines.split(/\r\n|\n/)
         //lines = lines.slice(1,-1);
         parse();
         plotdata(t,y,e);
      } else {
         document.getElementById("output").innerHTML = "<a href='"+url+"'>"+url+"</a><br> could not be loaded. Maybe a <a href='https://en.wikipedia.org/wiki/Cross-origin_resource_sharing'>CORS</a> issue. Try:<br>http://cdsarc.unistra.fr/ftp/J/A+A/552/A78/harps/hr3259_h.dat"
// https://raw.githubusercontent.com/plotly/datasets/master/spectral.csv
      }
   })
}

function handleFiles(files) {
   // Check for the various File API support.
   if (window.FileReader) {
      // FileReader are supported.
      getAsText(files[0]);
   } else {
      alert('FileReader are not supported in this browser.');
   }
}

function getAsText(fileToRead) {
   var reader = new FileReader();
   // Handle errors load
   reader.onload = loadHandler;
   reader.onerror = errorHandler;
   // Read file into memory as UTF-8
   reader.readAsText(fileToRead);
}

function loadHandler(event) {
   var csv = event.target.result;
   lines = csv.split(/\r\n|\n/);
   parse();
   plotdata(t,y);
}

function errorHandler(evt) {
   if(evt.target.error.name == "NotReadableError") {
      alert("Cannot read file !");
   }
}

