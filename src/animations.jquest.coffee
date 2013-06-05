$(window).load ->
  keyframes.setBrowserCode()
  keyframes.add [{
    name: "destinationRotate",
    "from": keyframes.browserCode+"transform:rotate(0deg)",
    "to": keyframes.browserCode+"transform:rotate(360deg)",
  }]