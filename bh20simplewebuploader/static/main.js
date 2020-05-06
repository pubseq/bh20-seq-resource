function fetchAPI(apiEndPoint) {
  fetch(scriptRoot + apiEndPoint)
    .then(response => {
      return response.json();
    })
    .then(data => {
      document.getElementById("json").textContent = JSON.stringify(data, undefined, 2);
      document.getElementById("results").classList.remove("invisible");
      document.getElementById("loader").classList.add("invisible");
    });
  document.getElementById("results").classList.add("invisible");
  document.getElementById("loader").classList.remove("invisible");

}

let search = () => {
  let m =  document.getElementById('search-input').value;
  fetchAPI(scriptRoot + "/api/getDetailsForSeq?seq=" + encodeURIComponent(m));
}

let fetchSEQBySpecimen = () => {
  fetchAPI("/api/getSEQCountbySpecimenSource");
}

let fetchSEQByLocation = () => {
  fetchAPI("/api/getSEQCountbyLocation");
}

let fetchSEQByTech = () => {
  fetchAPI("/api/getSEQCountbytech");
}

let fetchAllaccessions = () => {
  fetchAPI("/api/getAllaccessions");
};

/*
 * Make sure that only one of the manual metadata entry and metadata upload
 * form components is *actually* a child of the form element in the DOM.
 *
 * Because both make use of the "required" attribute, we can't get away with
 * just hiding the one we don't want the user to fill in. The hidden part will
 * still have possibly empty required fields and (some) browsers will
 * blocksubmission because of it. Moreover, the data (including file uploads)
 * from the hidden elements will still be sent to the server, which the user
 * may not expect.
 */

let uploadForm = document.getElementById('metadata_upload_form')
let uploadFormSpot = document.getElementById('metadata_upload_form_spot')
let fillForm = document.getElementById('metadata_fill_form')
let fillFormSpot = document.getElementById('metadata_fill_form_spot') 

function setUploadMode() {
  // Make the upload form the one in use.
  uploadFormSpot.appendChild(uploadForm)
  // Remove the upload form from the DOM so its required-ness does not block submission.
  fillFormSpot.removeChild(fillForm)
}

function setFillMode() {
  // Make the fillable form the one in use
  uploadFormSpot.removeChild(uploadForm)
  // Remove the fillable form from the DOM so its required-ness does not block submission.
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

/*
 * Machinery for variable-length lists of input items.
 */

// Start in mode appropriate to selected form item.
// It is important that we run this code when the page starts! The browser may
// have set the radio button to whatever the state was on last page load,
// instead of the default state, without raising an event, and we have to
// handle that.
setMode()

/**
 * Add another form field to the group this button is part of.
 */
function addField(e) {
  // Find our parent field-group div
  let fieldGroup = this.parentElement
  
  // Get its keypath
  let keypath = fieldGroup.dataset.keypath
  
  // Find its last field child
  let existingFields = fieldGroup.getElementsByClassName('field')
  let templateField = existingFields[existingFields.length - 1]
  
  // Get its number
  let fieldNumber = templateField.dataset.number
  
  // Duplicate it
  let newField = templateField.cloneNode(true)
  
  // Increment the number and use the keypath and number to set IDs and cross
  // references.
  // TODO: Heavily dependent on the form field HTML. Maybe we want custom
  // elements for the labeled controlsd that know how to be list items?
  fieldNumber++
  newField.dataset.number = fieldNumber
  let newID = keypath + '[' + fieldNumber + ']'
  let newControl = newField.getElementsByClassName('control')[0]
  newControl.id = newID
  newControl.setAttribute('name', newID)
  let newLabel = newField.getElementsByTagName('label')[0]
  newLabel.setAttribute('for', newID)
  
  // Find the minus button
  let minusButton = fieldGroup.getElementsByClassName('remove-field')[0]
  
  // Put new field as a child before the minus button
  fieldGroup.insertBefore(newField, minusButton)
  
  // Enable the minus button
  minusButton.classList.remove('invisible')
}

/**
 * Remove the last form field from the group button is part of.
 */
function removeField(e) {
  // Find our parent field-group div
  let fieldGroup = this.parentElement
  
  // Find its field children
  let existingFields = fieldGroup.getElementsByClassName('field')
  
  if (existingFields.length > 1) {
    // There is a last field we can safely remove.
    let lastField = existingFields[existingFields.length - 1]
    fieldGroup.removeChild(lastField)
  }
  
  if (existingFields.length <= 1) {
    // Collection auto-updates. Now there's only one element. Don't let the
    // user remove it. If they don't want it, they can leave it empty.
    this.classList.add('invisible')
  }
}

// Find all the add and remove field buttons and hook up the listeners.
for (let button of document.getElementsByClassName('add-field')) {
  button.addEventListener('click', addField)
}
for (let button of document.getElementsByClassName('remove-field')) {
  button.addEventListener('click', removeField)
}

