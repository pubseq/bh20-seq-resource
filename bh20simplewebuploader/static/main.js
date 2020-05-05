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

/**
 * Show form if checked
 */
let fillFormSpot = document.getElementById('metadata_fill_form_spot');
function displayForm() {
  if (document.getElementById('metadata_form').checked) {
    fillFormSpot.classList.remove("invisible");
    return;
  }
  fillFormSpot.classList.add("invisible");
}

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
  minusButton.classList.remove('hidden')
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
    this.classList.add('hidden')
  }
}

// Find all the add and remove field buttons and hook up the listeners.
for (let button of document.getElementsByClassName('add-field')) {
  button.addEventListener('click', addField)
}
for (let button of document.getElementsByClassName('remove-field')) {
  button.addEventListener('click', removeField)
}

