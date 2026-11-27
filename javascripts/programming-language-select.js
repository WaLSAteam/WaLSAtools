$(document).ready(function () {

    var tprolan = __md_get("__palette") || { index: e.findIndex((r) => matchMedia(r.getAttribute("data-md-color-media")).matches) };
    if (tprolan.index == 0) {
        var selectedplan = "python";
    } else if (tprolan.index == 1) {
        var selectedplan = "idl";
    }

    var urlstrg = "{{ nav_item.url }}";
    if (!urlstrg.match(selectedplan)) {
    document.getElementById("tabshtml").innerHTML +=  "<li class='md-tabs__item'><a href='{{ sitedomain+nav_item.url }}' class='{{ class }}'>{{ title }}</a></li>";
    }

    var urlstra = "{{ index.url }}";
    if (!urlstra.match(selectedplan)) {
        document.getElementById("navhtml1").innerHTML += "<div class='md-nav__link md-nav__link--index {{ class }}'><a href='{{ sitedomain+index.url }}'>{{ nav_item.title }}</a>{% if nav_item.children | length > 1 %}<label for='{{ path }}'><span class='md-nav__icon md-icon'></span></label>{% endif %}</div>";
    }

    var urlstrb = "{{ nav_item.url }}";
    if (!urlstrb.match(selectedplan)) {
    document.getElementById("navhtml2").innerHTML += "<a href='{{ sitedomain+nav_item.url }}' class='md-nav__link md-nav__link--active'>{{ nav_item.title }}</a>";
    }

    var urlstrc = "{{ nav_item.url }}";
    if (!urlstrc.match(selectedplan)) {
    document.getElementById("navhtml3").innerHTML += "<li class='{{ class }}'><a href='{{ sitedomain+nav_item.url }}' class='md-nav__link'>{{ nav_item.title }}</a></li>";
    }

});



