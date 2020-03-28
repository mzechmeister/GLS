var lines;

function handleFiles(files) {
   // Check for the various File API support.
   if (window.FileReader) {
      // FileReader are supported.
      getAsText(files[0]);
   } else {
      alert('FileReader are not supported in this browser.');
   }
}

function handleURL(url) {
   // Check for the various File API support.
   console.log(url)
   if (window.FileReader) {
      // FileReader are supported.
      getAsText(url);
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

