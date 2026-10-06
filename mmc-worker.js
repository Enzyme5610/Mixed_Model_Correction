// Reports webR downloads to the loading screen, then runs webR's worker.
// Copied next to webr-worker.js by the deploy workflow.
var bc = new BroadcastChannel("mmc-load");
function tell(state, u) { try { bc.postMessage({ state: state, url: String(u && u.url || u) }); } catch (e) {} }
var f0 = self.fetch;
self.fetch = function (u) {
  tell("start", u);
  return f0.apply(this, arguments).finally(function () { tell("done", u); });
};
var o0 = XMLHttpRequest.prototype.open;
XMLHttpRequest.prototype.open = function (m, u) {
  tell("start", u);
  this.addEventListener("loadend", function () { tell("done", u); });
  return o0.apply(this, arguments);
};
importScripts("webr-worker.js");
