$('.md-content').each(function () {
    var oldUrl = $(this).attr("href");
    if (tprolan.index == 0) {
        var newUrl = oldUrl.replace("/idl/", "/python/");
        $(this).attr("href", newUrl);
    } else if (tprolan.index == 1) {
        var newUrl = oldUrl.replace("/python/", "/idl/");
        $(this).attr("href", newUrl); 
    }
});