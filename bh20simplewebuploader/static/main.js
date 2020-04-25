let r = new XMLHttpRequest();
let test;
r.open("GET", scriptRoot + "/api/getAllaccessions", true);
r.onreadystatechange = function () {
  if (r.readyState != 4 || r.status != 200) return;
  test = r.responseText;
  console.log(JSON.parse(test));
};
r.send();
let uploadForm = document.getElementById('metadata_upload_form')
let uploadFormSpot = document.getElementById('metadata_upload_form_spot')
let fillForm = document.getElementById('metadata_fill_form')
let fillFormSpot = document.getElementById('metadata_fill_form_spot')

function setUploadMode() {
  // Make the upload form the one in use
  uploadFormSpot.appendChild(uploadForm)
  fillFormSpot.removeChild(fillForm)
}

function setFillMode() {
  // Make the fillable form the one in use
  uploadFormSpot.removeChild(uploadForm)
  fillFormSpot.appendChild(fillForm)
}

function setMode() {
  // Pick mode based on radio
  if (document.getElementById('metadata_upload').checked) {
    setUploadMode()
  } else {
    setFillMode()
  }
}

// Start in mode appropriate to selected form item
setMode()
